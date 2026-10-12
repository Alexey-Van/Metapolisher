#!/usr/bin/env python3
"""
partition_reads.py — разбивка ридов (ONT/HiFi/Illumina) по двум гаплотипам
фазированной сборки по разнице alignment score (AS) между выравниваниями
на hap1 и hap2.

Источник ридов (одно из двух):
  --reads FILE        один файл (ONT, HiFi или single-end Illumina)
  --r1 R1 --r2 R2     парные WGS-риды Illumina в двух файлах. Пара всегда
                      целиком попадает в одну категорию, порядок ридов в R1 и
                      R2 сохраняется (файлы остаются синхронными).

Категории (каждая в своём файле, чтобы доля неразделённых ридов была видна):
  hap1       уверенно лучше выровнялся на hap1 (AS1 - AS2 >= min_diff)
  hap2       уверенно лучше выровнялся на hap2 (AS2 - AS1 >= min_diff)
  shared     счёты почти равны, выравнивание однозначное (MAPQ >= min_mapq):
             рид из участка, где гаплотипы не различаются
  ambiguous  счёты почти равны и MAPQ низкий на обоих (повторы)
  unmapped   рид не выровнялся ни на один гаплотип

Для полишинга hap1 обычно берут: hap1 + shared (+ ambiguous по желанию).

Оценка точности (--stats): если истинный гаплотип закодирован в имени рида
(например, @chr1_PATERNAL:...:art_modern:8:nompi:62388/1), то после разбивки
печатается confusion-матрица (истина x предсказание) и метрики: сколько
ридов присвоено верно/неверно, accuracy среди присвоенных, доля присвоенных.
Какая метка соответствует hap1/hap2, задаётся --truth-hap1/--truth-hap2.

Вход: два BAM (риды vs hap1, риды vs hap2) с тегом AS. Считаются все
не-secondary записи рида (primary + supplementary), AS суммируется.
Для парных ридов суммируются оба мейта, поэтому BAM должны быть получены
из тех же R1/R2 (оба мейта в одном BAM).
"""

import argparse
import gzip
import sys
import zlib
from collections import Counter
from itertools import zip_longest

import pysam

CATS = ["hap1", "hap2", "shared", "ambiguous", "unmapped"]
HIST_CLIP = 500
TRUTHS = ["hap1", "hap2", "unknown"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Partition reads into hap1/hap2 by alignment score difference"
    )
    p.add_argument("--bam1", required=True, help="BAM: риды vs hap1")
    p.add_argument("--bam2", required=True, help="BAM: риды vs hap2")
    p.add_argument("--reads", help="Риды одним файлом (.fa/.fq, можно .gz)")
    p.add_argument("--r1", help="Парные риды: R1 (.fq/.fastq, можно .gz)")
    p.add_argument("--r2", help="Парные риды: R2 (.fq/.fastq, можно .gz)")
    p.add_argument("--out-prefix", required=True, help="Префикс выходных файлов")
    p.add_argument("--min-diff", type=int, default=1,
                   help="Минимальная абсолютная разница AS для присвоения "
                        "гаплотипа. Подбирать по --diff-hist (default: 10)")
    p.add_argument("--min-mapq", type=int, default=1,
                   help="MAPQ, ниже которого равный по счёту рид считается "
                        "ambiguous, а не shared (default: 1)")
    p.add_argument("--passes", type=int, default=1,
                   help="Число проходов по BAM/ридам для ограничения памяти: "
                        "в каждом проходе обрабатывается 1/N имён ридов "
                        "(default: 1)")
    p.add_argument("--threads", type=int, default=4,
                   help="Потоки для распаковки BAM (default: 4)")
    p.add_argument("--gzip", action="store_true", help="Сжимать выходные файлы")
    p.add_argument("--diff-hist", default=None,
                   help="Записать TSV-гистограмму разницы AS1-AS2 "
                        "(для выбора --min-diff)")
    p.add_argument("--stats", action="store_true",
                   help="Посчитать точность разбивки по меткам гаплотипа в "
                        "именах ридов; пишет также {out-prefix}.stats.tsv")
    p.add_argument("--truth-hap1", default="PATERNAL",
                   help="Подстрока имени рида, означающая, что рид из hap1 "
                        "(default: PATERNAL)")
    p.add_argument("--truth-hap2", default="MATERNAL",
                   help="Подстрока имени рида, означающая, что рид из hap2 "
                        "(default: MATERNAL)")
    args = p.parse_args()

    if args.stats and args.truth_hap1 == args.truth_hap2:
        p.error("--truth-hap1 и --truth-hap2 должны различаться")
    paired = args.r1 is not None or args.r2 is not None
    if paired and (args.r1 is None or args.r2 is None):
        p.error("--r1 и --r2 нужно указывать вместе")
    if paired and args.reads:
        p.error("используйте либо --reads, либо --r1/--r2, но не оба")
    if not paired and not args.reads:
        p.error("нужен источник ридов: --reads или --r1 и --r2")
    if args.passes < 1:
        p.error("--passes должен быть >= 1")
    return args


def norm_name(name):
    """Первое слово заголовка без суффикса /1 или /2."""
    name = name.split()[0]
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    return name


def in_pass(name, passes, k):
    return passes == 1 or zlib.crc32(name.encode()) % passes == k


def load_scores(bam_path, passes, k, threads):
    """name -> [сумма AS по не-secondary записям, максимальный MAPQ]."""
    scores = {}
    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary:
                continue
            name = norm_name(r.query_name)
            if not in_pass(name, passes, k):
                continue
            try:
                a = r.get_tag("AS")
            except KeyError:
                sys.exit(f"Ошибка: у записи {r.query_name} в {bam_path} нет "
                         f"тега AS. Нужен aligner, который его пишет "
                         f"(minimap2, bwa-mem2).")
            cur = scores.get(name)
            if cur is None:
                scores[name] = [a, r.mapping_quality]
            else:
                cur[0] += a
                cur[1] = max(cur[1], r.mapping_quality)
    return scores


def classify(s1, s2, min_diff, min_mapq):
    """Возвращает (категория, разница AS1-AS2 или None)."""
    if s1 is None and s2 is None:
        return "unmapped", None
    a1 = s1[0] if s1 else 0
    a2 = s2[0] if s2 else 0
    d = a1 - a2
    if abs(d) >= min_diff:
        return ("hap1" if d > 0 else "hap2"), d
    best_q = max(s[1] for s in (s1, s2) if s)
    return ("shared" if best_q >= min_mapq else "ambiguous"), d


def is_fastq(path):
    p = path[:-3] if path.endswith(".gz") else path
    return p.endswith((".fq", ".fastq"))


def open_in(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_fasta(fh):
    header, parts = None, []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(parts), None
            header, parts = line[1:], []
        else:
            parts.append(line)
    if header is not None:
        yield header, "".join(parts), None


def iter_fastq(fh):
    while True:
        h = fh.readline()
        if not h:
            break
        seq = fh.readline().rstrip("\n")
        fh.readline()
        qual = fh.readline().rstrip("\n")
        yield h[1:].rstrip("\n"), seq, qual


def iter_file(path, fh):
    return iter_fastq(fh) if is_fastq(path) else iter_fasta(fh)


def iter_units(args):
    """Выдаёт по одной «единице» на рид: [rec] или пару [rec1, rec2]."""
    if args.reads:
        with open_in(args.reads) as fh:
            for rec in iter_file(args.reads, fh):
                yield [rec]
        return
    with open_in(args.r1) as f1, open_in(args.r2) as f2:
        n = 0
        for rec1, rec2 in zip_longest(iter_file(args.r1, f1),
                                      iter_file(args.r2, f2)):
            n += 1
            if rec1 is None or rec2 is None:
                sys.exit(f"Ошибка: R1 и R2 разной длины (расхождение после "
                         f"{n - 1:,} пар).")
            if norm_name(rec1[0]) != norm_name(rec2[0]):
                sys.exit(f"Ошибка: имена R1/R2 не совпадают на паре {n:,}: "
                         f"{rec1[0].split()[0]} vs {rec2[0].split()[0]}. "
                         f"Файлы должны быть в одинаковом порядке.")
            yield [rec1, rec2]


def write_read(fh, header, seq, qual):
    if qual is None:
        fh.write(f">{header}\n{seq}\n")
    else:
        fh.write(f"@{header}\n{seq}\n+\n{qual}\n")


def out_paths(args):
    """Для каждой категории: список путей (1 файл или R1+R2)."""
    gz = ".gz" if args.gzip else ""

    def ext(path):
        return ("fq" if is_fastq(path) else "fa") + gz

    paths = {}
    for c in CATS:
        if args.reads:
            paths[c] = [f"{args.out_prefix}.{c}.{ext(args.reads)}"]
        else:
            paths[c] = [f"{args.out_prefix}.{c}.R1.{ext(args.r1)}",
                        f"{args.out_prefix}.{c}.R2.{ext(args.r2)}"]
    return paths


def truth_of(header, args):
    """Истинный гаплотип по имени рида: 'hap1', 'hap2' или 'unknown'."""
    name = header.split()[0]
    in1 = args.truth_hap1 in name
    in2 = args.truth_hap2 in name
    if in1 and not in2:
        return "hap1"
    if in2 and not in1:
        return "hap2"
    return "unknown"


def ratio(a, b):
    return a / b if b else None


def fmt_ratio(x):
    return "n/a" if x is None else f"{x:.4%}"


def report_stats(conf, args, unit_name):
    """Печатает confusion-матрицу и метрики, пишет {prefix}.stats.tsv."""
    total = {t: sum(conf[(t, c)] for c in CATS) for t in TRUTHS}
    known = total["hap1"] + total["hap2"]

    ok1, ok2 = conf[("hap1", "hap1")], conf[("hap2", "hap2")]
    bad1, bad2 = conf[("hap1", "hap2")], conf[("hap2", "hap1")]
    correct, wrong = ok1 + ok2, bad1 + bad2
    assigned = correct + wrong

    metrics = [
        (f"{'reads' if args.reads else 'pairs'}_with_known_truth", f"{known:,}"),
        ("assigned_correct", f"{correct:,}"),
        ("assigned_wrong", f"{wrong:,}"),
        ("unassigned_(shared+ambiguous+unmapped)", f"{known - assigned:,}"),
        ("accuracy_among_assigned", fmt_ratio(ratio(correct, assigned))),
        ("assignment_rate", fmt_ratio(ratio(assigned, known))),
        ("correct_rate_overall", fmt_ratio(ratio(correct, known))),
        ("wrong_rate_overall", fmt_ratio(ratio(wrong, known))),
        ("hap1_true_assigned_correctly", fmt_ratio(ratio(ok1, total["hap1"]))),
        ("hap2_true_assigned_correctly", fmt_ratio(ratio(ok2, total["hap2"]))),
        ("predicted_hap1_precision", fmt_ratio(ratio(ok1, ok1 + bad2))),
        ("predicted_hap2_precision", fmt_ratio(ratio(ok2, ok2 + bad1))),
    ]

    head = ["truth\\pred", "total"] + CATS
    rows = [head]
    for t in TRUTHS:
        if t == "unknown" and total[t] == 0:
            continue
        rows.append([t, f"{total[t]:,}"] + [f"{conf[(t, c)]:,}" for c in CATS])

    print(f"\nСтатистика точности (истина из имён; {unit_name}):",
          file=sys.stderr)
    for r in rows:
        print("  " + "".join(f"{x:>14}" for x in r), file=sys.stderr)
    print(file=sys.stderr)
    for k, v in metrics:
        print(f"  {k:<42}{v:>14}", file=sys.stderr)

    if known == 0:
        print(f"\nВНИМАНИЕ: в именах не найдено меток '{args.truth_hap1}' / "
              f"'{args.truth_hap2}', статистика пуста.", file=sys.stderr)
    elif assigned and wrong > correct:
        print(f"\nВНИМАНИЕ: неверных присвоений больше, чем верных. Возможно, "
              f"hap1 соответствует '{args.truth_hap2}', а не "
              f"'{args.truth_hap1}': попробуйте --truth-hap1 "
              f"{args.truth_hap2} --truth-hap2 {args.truth_hap1}.",
              file=sys.stderr)

    path = f"{args.out_prefix}.stats.tsv"
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(x.replace(",", "") for x in r) + "\n")
        fh.write("\n")
        for k, v in metrics:
            fh.write(f"{k}\t{v.replace(',', '')}\n")
    print(f"\n  Статистика записана в {path}", file=sys.stderr)


def main():
    args = parse_args()
    paired = args.reads is None
    opener = (lambda p: gzip.open(p, "wt")) if args.gzip else (lambda p: open(p, "w"))

    paths = out_paths(args)
    outs = {c: [opener(p) for p in ps] for c, ps in paths.items()}
    counts = Counter()
    hist = Counter()
    conf = Counter()

    try:
        for k in range(args.passes):
            tag = f"[проход {k + 1}/{args.passes}]"
            print(f"{tag} читаю {args.bam1}", file=sys.stderr)
            sc1 = load_scores(args.bam1, args.passes, k, args.threads)
            print(f"{tag} читаю {args.bam2}", file=sys.stderr)
            sc2 = load_scores(args.bam2, args.passes, k, args.threads)

            assignment = {}
            for name in sc1.keys() | sc2.keys():
                cat, d = classify(sc1.get(name), sc2.get(name),
                                  args.min_diff, args.min_mapq)
                assignment[name] = cat
                if d is not None:
                    hist[max(-HIST_CLIP, min(HIST_CLIP, d))] += 1
            del sc1, sc2

            print(f"{tag} записываю риды", file=sys.stderr)
            for unit in iter_units(args):
                name = norm_name(unit[0][0])
                if not in_pass(name, args.passes, k):
                    continue
                cat = assignment.get(name, "unmapped")
                if args.stats:
                    conf[(truth_of(unit[0][0], args), cat)] += 1
                for (header, seq, qual), fh in zip(unit, outs[cat]):
                    write_read(fh, header, seq, qual)
                counts[cat] += 1
            del assignment
    finally:
        for fhs in outs.values():
            for fh in fhs:
                fh.close()

    if args.diff_hist:
        with open(args.diff_hist, "w") as fh:
            fh.write("diff_AS1_minus_AS2\tn_reads\n")
            for d in sorted(hist):
                fh.write(f"{d}\t{hist[d]}\n")

    total = sum(counts.values())
    unit_name = "пар" if paired else "ридов"
    print(f"\nГотово ({unit_name}):", file=sys.stderr)
    for c in CATS:
        frac = counts[c] / total if total else 0
        print(f"  {c:<10} {counts[c]:>12,} ({frac:.1%})  "
              f"{' '.join(paths[c])}", file=sys.stderr)

    if args.stats:
        report_stats(conf, args, unit_name)


if __name__ == "__main__":
    main()