#!/usr/bin/env python3
import argparse
import gzip
import subprocess
from pathlib import Path


def sam_groups(bam):
    cmd = ['samtools', 'view', '-F', '2304', str(bam)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    current = None
    records = []

    try:
        for line in proc.stdout:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                continue

            qname = fields[0]
            if current is not None and qname != current:
                yield current, records
                records = []
            current = qname
            records.append(fields)

        if current is not None:
            yield current, records
    finally:
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f'samtools view failed for {bam} with exit code {rc}')


def score_group(records):
    mapped = []
    for f in records:
        flag = int(f[1])
        if flag & 4:
            continue
        mapq = int(f[4])
        tags = {}
        for field in f[11:]:
            p = field.split(':', 2)
            if len(p) == 3:
                tags[p[0]] = p[2]
        if 'AS' not in tags:
            continue
        try:
            as_score = int(tags['AS'])
            nm = int(tags.get('NM', '0'))
        except ValueError:
            continue
        mapped.append((as_score, mapq, nm))

    if not mapped:
        return None

    return (
        sum(x[0] for x in mapped),
        min(x[1] for x in mapped),
        sum(x[2] for x in mapped),
    )


def fastq_record(fields):
    flag = int(fields[1])
    seq = fields[9]
    qual = fields[10]

    if seq == '*':
        return None

    if flag & 16:
        complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
        seq = seq.translate(complement)[::-1]
        qual = qual[::-1]

    name = fields[0]
    return f'@{name}\n{seq}\n+\n{qual}\n'


def classify_and_write(hap1_bam, hap2_bam, delta, min_mapq, outdir, paired):
    it1 = sam_groups(hap1_bam)
    it2 = sam_groups(hap2_bam)
    a = next(it1, None)
    b = next(it2, None)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    classes = ('hap1', 'hap2', 'shared', 'unassigned')
    handles = {}

    for cls in classes:
        if paired:
            handles[(cls, 1)] = gzip.open(outdir / f'{cls}_R1.fastq.gz', 'wt')
            handles[(cls, 2)] = gzip.open(outdir / f'{cls}_R2.fastq.gz', 'wt')
        else:
            handles[cls] = gzip.open(outdir / f'{cls}.fastq.gz', 'wt')

    stats = {cls: 0 for cls in classes}

    try:
        while a is not None or b is not None:
            if b is None or (a is not None and a[0] < b[0]):
                name = a[0]
                records1 = a[1]
                records2 = []
                a = next(it1, None)
            elif a is None or b[0] < a[0]:
                name = b[0]
                records1 = []
                records2 = b[1]
                b = next(it2, None)
            else:
                name = a[0]
                records1 = a[1]
                records2 = b[1]
                a = next(it1, None)
                b = next(it2, None)

            s1 = score_group(records1)
            s2 = score_group(records2)

            if s1 is None and s2 is None:
                cls = 'unassigned'
                source = records1 or records2
            elif s1 is None:
                cls = 'hap2' if s2[1] >= min_mapq else 'unassigned'
                source = records2
            elif s2 is None:
                cls = 'hap1' if s1[1] >= min_mapq else 'unassigned'
                source = records1
            else:
                as1, mapq1, nm1 = s1
                as2, mapq2, nm2 = s2
                denom = max(abs(as1), abs(as2), 1)
                d = (as1 - as2) / denom

                if d >= delta and mapq1 >= min_mapq:
                    cls = 'hap1'
                    source = records1
                elif d <= -delta and mapq2 >= min_mapq:
                    cls = 'hap2'
                    source = records2
                elif abs(d) < delta and mapq1 >= min_mapq and mapq2 >= min_mapq:
                    cls = 'shared'
                    source = records1
                else:
                    cls = 'unassigned'
                    source = records1 or records2

            stats[cls] += 1

            if cls == 'unassigned':
                continue

            for f in source:
                flag = int(f[1])
                if flag & 2048 or flag & 256:
                    continue
                rec = fastq_record(f)
                if rec is None:
                    continue

                if paired:
                    if flag & 64:
                        handles[(cls, 1)].write(rec)
                    elif flag & 128:
                        handles[(cls, 2)].write(rec)
                else:
                    handles[cls].write(rec)

    finally:
        for h in handles.values():
            h.close()

    with open(outdir / 'summary.tsv', 'w') as out:
        out.write('class\treads_or_pairs\n')
        for cls in classes:
            out.write(f'{cls}\t{stats[cls]}\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hap1-bam', required=True)
    ap.add_argument('--hap2-bam', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--delta', type=float, default=0.02)
    ap.add_argument('--min-mapq', type=int, default=20)
    ap.add_argument('--paired', action='store_true')
    args = ap.parse_args()

    if not 0 <= args.delta < 1:
        raise ValueError('--delta must be >= 0 and < 1')

    classify_and_write(
        args.hap1_bam,
        args.hap2_bam,
        args.delta,
        args.min_mapq,
        args.outdir,
        args.paired,
    )


if __name__ == '__main__':
    main()
