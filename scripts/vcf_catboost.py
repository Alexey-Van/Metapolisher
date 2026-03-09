#!/usr/bin/env python3
"""
Variant ensemble ML pipeline with TRUE ground truth (SV-aware)

- Merge multiple VCFs into a single variant table
- Annotate with GFF3 / BED (RepeatMasker, Liftoff, low complexity, Flagger)
- Build TRUE_VARIANT label from truth VCF (±1bp)
- Split by SVTYPE (SNV / INDEL)
- Train CatBoost models independently
- Save models
- Write SNV / INDEL / ALL filtered VCFs
"""

import argparse
import pysam
from pathlib import Path
import pandas as pd
import numpy as np
from intervaltree import IntervalTree
from cyvcf2 import VCF, Writer
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return 0.0


def normalize_info_value(key, value):
    if value is None:
        return 0
    if key.upper() == "SVLEN":
        if isinstance(value, (list, tuple)):
            return abs(safe_float(value[0]))
        return abs(safe_float(value))
    if isinstance(value, (list, tuple)):
        try:
            return float(np.mean(value))
        except Exception:
            return 0
    return safe_float(value)


def normalize_svtype(x):
    if x in ("SNV", "SNP"):
        return "SNV"
    return "INDEL"


# ------------------------------------------------------------
# CatBoost
# ------------------------------------------------------------

def train_catboost(df, label="TRUE_VARIANT"):
    drop_cols = {
        label, "CHROM", "POS", "END",
        "REF", "ALT", "SVTYPE", "SVTYPE_NORM"
    }

    X = df.drop(columns=[c for c in drop_cols if c in df.columns])
    y = df[label]

    if y.nunique() < 2:
        raise RuntimeError("Only one class present — training impossible")

    cat_cols = X.select_dtypes(include="object").columns.tolist()

    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, test_size=0.8, stratify=y, random_state=42
    )

    model = CatBoostClassifier(
        iterations=300,
        depth=4,
        learning_rate=0.05,
        loss_function="Logloss", 
        auto_class_weights="Balanced",
        l2_leaf_reg=5,
        eval_metric="AUC",
        verbose=100,
    )   


    model.fit(X_tr, y_tr, cat_features=cat_cols, eval_set=(X_te, y_te))

    # после model.fit(...)
    print("\n=== Feature Importance (PredictionValuesChange) ===")
    fi = model.get_feature_importance(type="PredictionValuesChange")
    for name, val in sorted(zip(X.columns, fi), key=lambda x: -x[1])[:30]:
        print(f"{name:40s} {val:.5f}")

    auc = roc_auc_score(y_te, model.predict_proba(X_te)[:, 1])
    print(f"AUC = {auc:.4f}")
    

    return model


def predict_block(df, model):
    drop_cols = {
        "TRUE_VARIANT", "CHROM", "POS", "END",
        "REF", "ALT", "SVTYPE", "SVTYPE_NORM"
    }
    X = df.drop(columns=[c for c in drop_cols if c in df.columns])
    return model.predict_proba(X)[:, 1]


# ------------------------------------------------------------
# FORMAT parsing
# ------------------------------------------------------------

def parse_format_fields(v, prefix):
    out = {f"{prefix}_QUAL": safe_float(v.QUAL)}

    if not v.FORMAT or not v.genotypes:
        return out

    for fmt in v.FORMAT:
        try:
            arr = v.format(fmt)
        except KeyError:
            continue
        if arr is None or len(arr) == 0:
            continue

        val = arr[0]
        key = f"{prefix}_FORMAT_{fmt}"

        if fmt == "GT" and isinstance(val, (list, np.ndarray)):
            out[f"{prefix}_GT"] = int(val[0] != val[1])
        elif isinstance(val, (int, float, np.integer, np.floating)):
            out[key] = float(val)
        elif isinstance(val, (list, np.ndarray)):
            nums = [x for x in val if isinstance(x, (int, float))]
            out[key] = float(np.mean(nums)) if nums else 0.0
        else:
            out[key] = str(val)

    return out


# ------------------------------------------------------------
# VCF loading
# ------------------------------------------------------------

def load_vcf_as_df(path):
    prefix = Path(path).stem
    rows = []

    for v in VCF(path):
        info = dict(v.INFO)

        row = {
            "CHROM": v.CHROM,
            "POS": v.POS,
            "END": max(v.end, v.POS),
            "REF": v.REF,
            "ALT": str(v.ALT[0]) if v.ALT else None,
            "SVTYPE": info.get("SVTYPE", "SNV"),
            "SVLEN": normalize_info_value("SVLEN", info.get("SVLEN")),
            f"{prefix}_present": 1,
        }

        for k, val in info.items():
            row[f"{prefix}_INFO_{k}"] = normalize_info_value(k, val)

        row.update(parse_format_fields(v, prefix))
        rows.append(row)

    return pd.DataFrame(rows)

def load_merfin_support(merfin_vcf):
    """
    Загружает Merfin PASS VCF и возвращает множество поддержанных позиций.
    Ключ: (CHROM, POS)
    """
    supported = set()
    for v in VCF(merfin_vcf):
        supported.add((v.CHROM, v.POS))
    return supported


def annotate_merfin_support(df, merfin_vcf):
    """
    Добавляет колонку MERFIN_SUPPORTED в общий DataFrame.
    """
    supported = load_merfin_support(merfin_vcf)

    df["MERFIN_SUPPORTED"] = df.apply(
        lambda r: int((r.CHROM, r.POS) in supported),
        axis=1
    )
    return df

def load_merfin_scores(bed_path):
    """
    Загружает Merfin kmer_scores.bed в IntervalTree по хромосомам.
    """
    trees = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, score = line.rstrip().split("\t")
            start, end = int(start), int(end)
            score = float(score)
            trees.setdefault(chrom, IntervalTree()).addi(start, end, score)
    return trees


def annotate_merfin_scores(df, bed_path):
    """
    Добавляет колонку MERFIN_SCORE (или 0.0, если нет попадания).
    """
    trees = load_merfin_scores(bed_path)
    scores = []

    for _, r in df.iterrows():
        hits = trees.get(r.CHROM, IntervalTree()).overlap(r.POS, r.END)
        if hits:
            scores.append(next(iter(hits)).data)
        else:
            scores.append(0.0)

    df["MERFIN_SCORE"] = scores
    return df


def outer_merge_vcfs(paths):
    dfs = [load_vcf_as_df(p) for p in paths]
    df = dfs[0]
    for d in dfs[1:]:
        df = df.merge(
            d,
            on=["CHROM", "POS", "END", "REF", "ALT", "SVTYPE", "SVLEN"],
            how="outer",
        )
    return df.fillna(0)

def collect_contigs(vcfs):
    contigs = set()
    for path in vcfs:
        for v in VCF(path):
            contigs.add(v.CHROM)
    return sorted(contigs)

def fixed_header(template_vcf, all_vcfs):
    vf = pysam.VariantFile(template_vcf)
    header = vf.header.copy()

    # --- contigs ---
    contigs = set()
    for vcf in all_vcfs:
        for rec in pysam.VariantFile(vcf).fetch():
            contigs.add(rec.chrom)

    for c in sorted(contigs):
        if c not in header.contigs:
            header.contigs.add(c)

    # --- INFO ---
    if "ML_PROB" not in header.info:
        header.info.add(
            "ML_PROB",
            number=1,
            type="Float",
            description="CatBoost predicted probability"
        )

    return header





# ------------------------------------------------------------
# Annotation loaders
# ------------------------------------------------------------

def load_bed(path, value_col=None):
    trees = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                if len(parts) < 3:
                    continue

            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            val = parts[value_col] if value_col is not None and len(parts) > value_col else 1
            trees.setdefault(chrom, IntervalTree()).addi(start, end, val)
    return trees



def load_gff3(path, attr_key=None):
    trees = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            if start == end:
                end += 1

            value = parts[2]
            if attr_key:
                attrs = dict(
                    x.split("=", 1) for x in parts[8].split(";") if "=" in x
                )
                raw = attrs.get(attr_key)
                value = raw.split(" ", 1)[0] if raw else None

            trees.setdefault(chrom, IntervalTree()).addi(start, end, value)

    return trees


# ------------------------------------------------------------
# Annotation application
# ------------------------------------------------------------

def annotate_interval(df, trees, colname):
    present, value = [], []

    for _, r in df.iterrows():
        hits = trees.get(r.CHROM, IntervalTree()).overlap(r.POS, r.END + 1)
        if hits:
            h = next(iter(hits))
            present.append(1)
            value.append(h.data)
        else:
            present.append(0)
            value.append("NA")

    df[f"{colname}_present"] = present
    df[colname] = value


# ------------------------------------------------------------
# Truth
# ------------------------------------------------------------

def load_truth(vcf):
    truth = {}
    for v in VCF(vcf):
        truth.setdefault(v.CHROM, set()).add(v.POS)
    return truth


def is_true(r, truth):
    return int(any(abs(r.POS - p) <= 1 for p in truth.get(r.CHROM, [])))


# ------------------------------------------------------------
# VCF writing
# ------------------------------------------------------------

def write_vcf(df, vcfs, out_vcf, threshold):
    """
    - keep only variants with ML_PROB >= threshold
    - if multiple rows share same (CHROM, POS): keep max ML_PROB
    - INFO:
        ML_PROB
        SOURCE=tool1,tool2,...
    """

    # ----------------------------
    # 1. filter by threshold
    # ----------------------------
    df_filt = df[df.ML_PROB >= threshold].copy()
    if df_filt.empty:
        print(f"[VCF] no variants pass threshold for {out_vcf}")
        return

    # ----------------------------
    # 2. deduplicate by CHROM, POS
    # ----------------------------
    df_filt.sort_values("ML_PROB", ascending=False, inplace=True)
    df_best = df_filt.drop_duplicates(subset=["CHROM", "POS"], keep="first")

    # ----------------------------
    # 3. prepare header
    # ----------------------------
    header = fixed_header(vcfs[0], vcfs)

    if "SOURCE" not in header.info:
        header.info.add(
            "SOURCE",
            number=".",
            type="String",
            description="Source tools contributing to this variant"
        )

    out = pysam.VariantFile(out_vcf, "w", header=header)

    # ----------------------------
    # 4. index variants by key
    # ----------------------------
    best = {}
    for _, r in df_best.iterrows():
        key = (r.CHROM, r.POS)
        best[key] = r

    # ----------------------------
    # 5. write records
    # ----------------------------
    written = 0

    for vcf in vcfs:
        vf = pysam.VariantFile(vcf)
        for v in vf.fetch():
            key = (v.chrom, v.pos)
            if key not in best:
                continue

            r = best[key]

            # ---- collect sources ----
            sources = []
            for c in df.columns:
                if c.endswith("_present") and r.get(c, 0) == 1:
                    sources.append(c.replace("_present", ""))

            nv = out.new_record(
                contig=v.chrom,
                start=v.start,
                stop=v.stop,
                alleles=v.alleles
            )

            nv.info["ML_PROB"] = float(r.ML_PROB)
            if sources:
                nv.info["SOURCE"] = ",".join(sorted(sources))

            out.write(nv)
            written += 1

    out.close()
    print(f"[VCF] written {written} variants → {out_vcf}")




# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main(args):
    print("[1] Loading and merging VCFs")
    df = outer_merge_vcfs(args.vcfs)
    df = annotate_merfin_support(df, "merfin.pass.vcf")
    df = annotate_merfin_scores(df, "kmer_scores.bed")


    print("[2] Loading annotations")
    repeat = load_gff3(args.repeat_gff, attr_key="Target")
    liftoff = load_gff3(args.liftoff_gff)
    low_complex = load_bed(args.low_complex)
    flagger = load_bed(args.flagger, value_col=3)

    print("[3] Annotating variants")
    annotate_interval(df, repeat, "repeat_target")
    annotate_interval(df, liftoff, "liftoff_feature")
    annotate_interval(df, low_complex, "low_complexity")
    annotate_interval(df, flagger, "flagger_state")

    print("[4] Building TRUE_VARIANT labels")
    df["SVTYPE_NORM"] = df["SVTYPE"].apply(normalize_svtype)
    truth = load_truth(args.truth_vcf)
    df["TRUE_VARIANT"] = df.apply(lambda r: is_true(r, truth), axis=1)

    print("[5] Splitting SNV / INDEL")
    df_snv = df[df.SVTYPE_NORM == "SNV"].copy()
    df_indel = df[df.SVTYPE_NORM == "INDEL"].copy()

    print(f"SNV   : {len(df_snv):,}")
    print(f"INDEL : {len(df_indel):,}")

    print("[6] Training models")
    model_snv = train_catboost(df_snv)
    model_snv.save_model("catboost_snv.cbm")

    model_indel = train_catboost(df_indel)
    model_indel.save_model("catboost_indel.cbm")

    print("[7] Predicting ML_PROB")
    df["ML_PROB"] = 0.0
    df.loc[df_snv.index, "ML_PROB"] = predict_block(df_snv, model_snv)
    df.loc[df_indel.index, "ML_PROB"] = predict_block(df_indel, model_indel)

    df_snv["ML_PROB"] = df.loc[df_snv.index, "ML_PROB"]
    df_indel["ML_PROB"] = df.loc[df_indel.index, "ML_PROB"]

    df.to_csv(args.out_table, sep="\t", index=False)


    print("[8] Writing VCFs")
    write_vcf(df_snv, args.vcfs, "snv.filtered.vcf", args.threshold)
    write_vcf(df_indel, args.vcfs, "indel.filtered.vcf", args.threshold)
    write_vcf(df, args.vcfs, args.out_vcf, args.threshold)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--vcfs", nargs="+", required=True)
    p.add_argument("--truth_vcf", required=True)
    p.add_argument("--merfin_pass_vcf", required=True)
    p.add_argument("--merfin_scores_bed", required=False)
    p.add_argument("--repeat_gff", required=True)
    p.add_argument("--liftoff_gff", required=True)
    p.add_argument("--low_complex", required=True)
    p.add_argument("--flagger", required=True)
    p.add_argument("--out_vcf", required=True)
    p.add_argument("--out_table", default="variant_features.tsv")
    p.add_argument("--threshold", type=float, default=0.9)
    main(p.parse_args())
