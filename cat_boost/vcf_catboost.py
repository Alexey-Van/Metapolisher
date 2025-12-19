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

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return 0.0

def safe_int(x):
    try:
        return int(x)
    except Exception:
        return 0

def normalize_info_value(key, value):
    if value is None:
        return 0
    if key.upper() == "COVERAGE":
        if isinstance(value, (list, tuple)) and value:
            return float(sum(value)) / len(value)
        return safe_float(value)
    if key.upper() == "SVLEN":
        if isinstance(value, (list, tuple)):
            return abs(safe_float(value[0]))
        return abs(safe_float(value))
    if isinstance(value, (list, tuple)):
        try:
            return float(sum(value)) / len(value)
        except Exception:
            return 0
    return safe_float(value) if isinstance(value, (int, float)) else value

# ------------------------------------------------------------
# FORMAT parsing
# ------------------------------------------------------------

def parse_format_fields(v, prefix):
    out = {}

    out[f"{prefix}_QUAL"] = float(v.QUAL) if v.QUAL is not None else 0.0

    if not v.FORMAT or not v.genotypes:
        return out

    sample_idx = 0 

    for fmt in v.FORMAT:
        try:
            arr = v.format(fmt)
        except KeyError:
            continue

        if arr is None or len(arr) == 0:
            continue

        val = arr[sample_idx]

        # ---------- GT ----------
        if fmt == "GT":
            if isinstance(val, (list, np.ndarray)) and len(val) >= 2:
                out[f"{prefix}_GT"] = int(val[0] != val[1])
            else:
                out[f"{prefix}_GT"] = 0
            continue

        # ---------- scalar numeric ----------
        if isinstance(val, (int, float, np.integer, np.floating)):
            out[f"{prefix}_FORMAT_{fmt}"] = float(val)
            continue

        # ---------- numeric array (AD, PL) ----------
        if isinstance(val, (list, np.ndarray)):
            nums = [x for x in val if isinstance(x, (int, float, np.integer, np.floating))]
            out[f"{prefix}_FORMAT_{fmt}"] = float(np.mean(nums)) if nums else 0.0
            continue

        # ---------- categorical string (FT=PASS etc.) ----------
        if isinstance(val, (str, np.str_)):
            out[f"{prefix}_FORMAT_{fmt}"] = str(val)
            continue

        # ---------- fallback ----------
        out[f"{prefix}_FORMAT_{fmt}"] = "NA"

    return out




# ------------------------------------------------------------
# VCF loading
# ------------------------------------------------------------

def load_vcf_as_df(vcf_path, prefix):
    records = []
    vcf = VCF(str(vcf_path))
    is_deepvariant = "deepvariant" in prefix.lower()

    for v in vcf:
        info = dict(v.INFO)
        
        pos = v.POS - 1 if is_deepvariant else v.POS
        end = v.end - 1 if (is_deepvariant and v.end) else v.end
        end = max(end, pos)

        row = {
            "CHROM": v.CHROM,
            "POS": v.POS,
            "END": max(v.end, v.POS),
            "REF": v.REF,
            "ALT": str(v.ALT[0]) if v.ALT else None,
            "SVTYPE": info.get("SVTYPE", "SNV"),
            "SVLEN": normalize_info_value("SVLEN", info.get("SVLEN", 0)),
            f"{prefix}_present": 1,
            f"{prefix}_QUAL": safe_float(v.QUAL),
        }

        # INFO
        for k, val in info.items():
            row[f"{prefix}_INFO_{k}"] = normalize_info_value(k, val)

        # FORMAT
        row.update(parse_format_fields(v, prefix))

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
    present = []
    value = []

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
# Truth VCF matching (±1bp)
# ------------------------------------------------------------

def load_truth_variants(truth_vcf):
    truth = {}
    for v in VCF(truth_vcf):
        truth.setdefault(v.CHROM, set()).add(v.POS)
    return truth

def is_true_variant(r, truth_dict):
    for tpos in truth_dict.get(r.CHROM, []):
        if abs(r.POS - tpos) <= 1:
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

    print("[4] Building TRUE_VARIANT labels")
    truth = load_truth_variants(args.truth_vcf)
    df["TRUE_VARIANT"] = df.apply(lambda r: is_true_variant(r, truth), axis=1)

    df.to_csv(args.out_table, sep="\t", index=False)

    print("[5] Training CatBoost")
    y = df.TRUE_VARIANT
    X = df.drop(columns=["TRUE_VARIANT"])
    cat_cols = X.select_dtypes(include="object").columns.tolist()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=0.8, random_state=42
    )

    model = CatBoostClassifier(
        iterations=500,
        depth=10,
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

    print("AUC:", roc_auc_score(y_test, model.predict_proba(X_test)[:, 1]))

    print("[6] Writing filtered VCF")
    base_vcf = VCF(args.vcfs[0])
    base_vcf.add_info_to_header({
        "ID": "ML_PROB",
        "Type": "Float",
        "Number": "1",
        "Description": "CatBoost TRUE_VARIANT probability",
    })

    writer = Writer(args.out_vcf, base_vcf)
    df["ML_PROB"] = model.predict_proba(X)[:, 1]

    prob_map = {
        (r.CHROM, r.POS, r.END, r.REF, r.ALT): r.ML_PROB
        for _, r in df.iterrows()
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
    main(p.parse_args())
