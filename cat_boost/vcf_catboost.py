#!/usr/bin/env python3
"""
Variant ensemble ML pipeline with TRUE ground truth (SV-aware).

- Merge multiple VCFs into a single variant table
- Annotate with GFF3 / BED (RepeatMasker, Liftoff, low complexity, Flagger)
- Build TRUE_VARIANT label from truth VCF (±1bp tolerance)
- Train CatBoostClassifier (binary classification)
- Output filtered VCF with ML probability
"""

import argparse
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

def load_vcf_as_df(vcf_path, prefix):
    records = []
    vcf = VCF(str(vcf_path))
    for v in vcf:
        info = dict(v.INFO)
        row = {
            "CHROM": v.CHROM,
            "POS": v.POS,
            "END": max(v.end, v.POS),  # ensure END >= POS
            "REF": v.REF,
            "ALT": str(v.ALT[0]) if v.ALT else None,
            "SVTYPE": info.get("SVTYPE", "SNV"),
            "SVLEN": info.get("SVLEN", 0),
            f"{prefix}_present": 1,
        }
        for k, val in info.items():
            row[f"{prefix}_INFO_{k}"] = normalize_info_value(k, val)
        records.append(row)
    return pd.DataFrame(records)

def outer_merge_vcfs(vcf_paths):
    dfs = []
    for p in vcf_paths:
        name = Path(p).stem
        dfs.append(load_vcf_as_df(p, name))
    df = dfs[0]
    for other in dfs[1:]:
        df = df.merge(
            other,
            on=["CHROM", "POS", "END", "REF", "ALT", "SVTYPE", "SVLEN"],
            how="outer",
        )
    df.fillna(0, inplace=True)
    return df

def normalize_info_value(key, value):
    if value is None:
        return None
    if key.upper() == "COVERAGE":
        if isinstance(value, (list, tuple)) and len(value) > 0:
            return float(sum(value))/len(value)
        try:
            return float(value)
        except Exception:
            return 0
    if key.upper() == "SVLEN":
        if isinstance(value, (list, tuple)):
            return abs(float(value[0]))
        return abs(float(value))
    if isinstance(value, (list, tuple)):
        try:
            return float(sum(value))/len(value)
        except Exception:
            return str(value)
    return value

# ------------------------------------------------------------
# Annotation loaders
# ------------------------------------------------------------

def load_bed(path, value_col=None):
    trees = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            val = parts[value_col] if value_col is not None and len(parts) > value_col else 1
            trees.setdefault(chrom, IntervalTree()).addi(start, end, val)
    return trees

def load_gff3(path, attr_key=None, col_idx=None):
    trees = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])

            if start > end:
                continue
            if start == end:
                end = start + 1  # fix zero-length intervals

            value = feature

            if attr_key:
                attrs = dict(
                    x.split("=", 1) for x in parts[8].split(";") if "=" in x
                )
                raw = attrs.get(attr_key)
                if raw:
                    # take everything before first space
                    value = raw.split(" ", 1)[0]
                else:
                    value = None

            if col_idx is not None and len(parts) > col_idx:
                value = parts[col_idx]

            trees.setdefault(chrom, IntervalTree()).addi(start, end, value)

    return trees



# ------------------------------------------------------------
# Annotation application
# ------------------------------------------------------------

def annotate_interval(df, trees, colname):
    values = []
    for _, r in df.iterrows():
        hits = trees.get(r.CHROM, IntervalTree()).overlap(r.POS, r.END+1)  # POS == END included
        values.append(next(iter(hits)).data if hits else None)
    df[colname] = values

# ------------------------------------------------------------
# Truth VCF matching (±1bp)
# ------------------------------------------------------------

def load_truth_variants(truth_vcf):
    truth = {}
    vcf = VCF(truth_vcf)
    for v in vcf:
        chrom, pos, end = v.CHROM, v.POS, max(v.end, v.POS)
        alt = str(v.ALT[0]) if v.ALT else None
        key = (chrom, pos, end, v.REF, alt)
        # store in dict for ±1 lookup
        truth.setdefault(chrom, set()).add((pos, end, v.REF, alt))
    return truth

def is_true_variant(r, truth_dict):
    chrom = r.CHROM
    if chrom not in truth_dict:
        return 0
    for t in truth_dict[chrom]:
        t_start, t_end, t_ref, t_alt = t
        if abs(r.POS - t_start) <= 1:
            return 1
    return 0

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main(args):
    print("[1] Loading and merging VCFs")
    df = outer_merge_vcfs(args.vcfs)
    print(f"Merged variants: {len(df):,}")

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

    # fill missing values to avoid CatBoost errors
    for col in df.columns:
        if df[col].dtype in [float, int]:
            df[col] = df[col].fillna(0)
        else:
            df[col] = df[col].fillna("NA")

    print("[4] Building TRUE_VARIANT labels")
    truth_dict = load_truth_variants(args.truth_vcf)
    df["TRUE_VARIANT"] = df.apply(lambda r: is_true_variant(r, truth_dict), axis=1)

    print("Saving feature table to TSV for inspection")
    df.to_csv(args.out_table, sep="\t", index=False)

    print("[5] Training CatBoost")
    y = df.TRUE_VARIANT
    X = df.drop(columns=["TRUE_VARIANT"])
    cat_cols = X.select_dtypes(include=["object"]).columns.tolist()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )

    model = CatBoostClassifier(
        iterations=500,
        depth=8,
        learning_rate=0.05,
        loss_function="Logloss",
        eval_metric="AUC",
        verbose=100
    )

    model.fit(
        X_train, y_train,
        cat_features=cat_cols,
        eval_set=(X_test, y_test)
    )

    preds = model.predict_proba(X_test)[:, 1]
    print("AUC:", roc_auc_score(y_test, preds))

    print("[6] Writing filtered VCF")
    base_vcf = VCF(args.vcfs[0])
    base_vcf.add_info_to_header({
        "ID": "ML_PROB",
        "Description": "CatBoost probability TRUE_VARIANT",
        "Type": "Float",
        "Number": "1",
    })
    writer = Writer(args.out_vcf, base_vcf)

    prob_map = {
        (r.CHROM, r.POS, r.END, r.REF, r.ALT): model.predict_proba(X.loc[[i]])[0][1]
        for i, r in df.iterrows()
    }

    for v in base_vcf:
        key = (v.CHROM, v.POS, max(v.end, v.POS), v.REF, str(v.ALT[0]) if v.ALT else None)
        if key in prob_map and prob_map[key] >= args.threshold:
            v.INFO["ML_PROB"] = round(prob_map[key], 4)
            writer.write_record(v)

    writer.close()
    print("Done")

# ------------------------------------------------------------

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--vcfs", nargs="+", required=True)
    p.add_argument("--truth_vcf", required=True)
    p.add_argument("--repeat_gff", required=True)
    p.add_argument("--liftoff_gff", required=True)
    p.add_argument("--low_complex", required=True)
    p.add_argument("--flagger", required=True)
    p.add_argument("--out_vcf", required=True)
    p.add_argument("--out_table", default="variant_features.tsv")
    p.add_argument("--threshold", type=float, default=0.9)
    args = p.parse_args()
    main(args)
