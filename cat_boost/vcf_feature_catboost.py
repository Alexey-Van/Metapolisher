#!/usr/bin/env python3
"""
vcf_feature_catboost.py

Usage:
  python vcf_feature_catboost.py --ref ref.fa --gff3 ann.gff3 --vcf sample1.vcf sample2.vcf \
    --repeat repeatmask.bed --flagger flagger.bed --bed genes.bed --out features.tsv --model out_model.cbm

Outputs:
  - features.tsv : table with features + TRUE_VCF label (if present in INFO)
  - out_model.cbm : trained CatBoost model
  - preds.csv : features + predicted probability
"""

import argparse
from pathlib import Path
import sys

import pandas as pd
import numpy as np

from pyfaidx import Fasta
import gffutils
from cyvcf2 import VCF
import pybedtools

from catboost import CatBoostClassifier, Pool
from sklearn.model_selection import GroupKFold
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def get_kmer_from_fa(fasta, chrom, pos0, k=11):
    """Return k-mer centered at pos0 (0-based). If out of bounds, return N..N"""
    half = k // 2
    start = pos0 - half
    end = pos0 + half + 1
    try:
        # pyfaidx slicing: fasta[chrom][start:end] returns str
        seq = fasta[chrom][max(0, start):max(0, start) + (end - max(0, start))].seq
        # if start < 0 or end > chrom_len => pad
        if start < 0:
            seq = "N" * (-start) + seq
        chrom_len = len(fasta[chrom])
        if end > chrom_len:
            seq = seq + "N" * (end - chrom_len)
        return seq.upper()
    except KeyError:
        return "N" * k

def gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return np.nan
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq)

def load_gff_db(gff3_path):
    db_path = str(gff3_path) + ".db"
    if not Path(db_path).exists():
        print(f"[INFO] Creating gffutils DB at {db_path} (this may take a while)...")
        gffutils.create_db(str(gff3_path), db_path, merge_strategy="merge", force=True)
    return gffutils.FeatureDB(db_path)

def region_and_gene_info(db, chrom, pos1):
    """pos1: 1-based position"""
    try:
        hits = list(db.region(region=(chrom, pos1, pos1)))
    except Exception:
        return "intergenic", None, None
    if not hits:
        return "intergenic", None, None
    types = {f.featuretype for f in hits}
    # prioritize
    if "CDS" in types:
        region = "CDS"
    elif "exon" in types:
        region = "exon"
    elif "intron" in types:
        region = "intron"
    elif "five_prime_UTR" in types or "three_prime_UTR" in types:
        region = "UTR"
    else:
        region = list(types)[0]
    gene_id = None
    gene_biotype = None
    for f in hits:
        if f.featuretype == "gene":
            gene_id = getattr(f, "id", None)
            gene_biotype = f.attributes.get("biotype", [None])[0] if hasattr(f, "attributes") else None
            break
    return region, gene_id, gene_biotype

def bed_any_hit(bedtool, chrom, pos1):
    """Return first hit.name or None. pos1 is 1-based."""
    if bedtool is None:
        return None
    # pybedtools uses 0-based half-open intervals; convert to 0-based start, end=pos
    interval = pybedtools.Interval(chrom, pos1 - 1, pos1)
    try:
        if bedtool.any_hits(interval):
            # find returns iterator of matching intervals; take first
            for hit in bedtool.intersect(pybedtools.BedTool([interval]), wao=True, u=True):
                # hit may be Interval; name at .name or 3rd/4th field; be defensive
                return getattr(hit, "name", None) or (hit.fields[3] if len(hit.fields) > 3 else None)
        return None
    except Exception:
        # fallback: try intersect with the whole bedtool
        try:
            hits = bedtool.intersect(pybedtools.BedTool([interval]))
            for h in hits:
                return getattr(h, "name", None) or (h.fields[3] if len(h.fields) > 3 else None)
            return None
        except Exception:
            return None

# ------------------------------------------------------------
# Feature extraction per VCF
# ------------------------------------------------------------
def extract_from_vcf(vcf_path, fasta, gff_db, bed_files=None, k=11, gc_window=101):
    """
    vcf_path: path to vcf
    fasta: pyfaidx Fasta object
    gff_db: gffutils DB
    bed_files: dict with keys: 'genes','repeat','flagger' mapped to pybedtools.BedTool or None
    """
    rows = []
    vcf = VCF(str(vcf_path))
    sample_name = Path(vcf_path).stem

    for v in vcf:
        chrom = v.CHROM
        pos1 = v.POS  # 1-based
        pos0 = pos1 - 1
        ref = v.REF
        alt = ",".join(v.ALT) if v.ALT is not None else ""
        var_type = v.var_type  # SNV, INDEL, etc.
        qual = v.QUAL if v.QUAL is not None else np.nan
        flt = ";".join(v.FILTER or []) if v.FILTER is not None else None

        # sequence context
        kmer = get_kmer_from_fa(fasta, chrom, pos0, k=k)
        win_start = max(0, pos0 - gc_window // 2)
        win_end = pos0 + gc_window // 2 + 1
        try:
            win_seq = fasta[chrom][win_start:win_end].seq
        except KeyError:
            win_seq = ""
        gc = gc_content(win_seq)

        # region/gene info from GFF3
        region, gene_id, gene_biotype = region_and_gene_info(gff_db, chrom, pos1)

        # bed annotations (genes/other)
        bed_annot = bed_any_hit(bed_files.get("genes") if bed_files else None, chrom, pos1)
        repeat_class = bed_any_hit(bed_files.get("repeat") if bed_files else None, chrom, pos1)
        flagger_state = bed_any_hit(bed_files.get("flagger") if bed_files else None, chrom, pos1)

        # collect INFO fields (flatten simple scalars; arrays -> join)
        info_dict = {}
        try:
            for kinfo, val in v.INFO.items():
                # cyvcf2 returns arrays or scalars; ensure stringifiable
                if isinstance(val, (list, tuple)):
                    info_dict[f"INFO_{kinfo}"] = ";".join(map(str, val))
                else:
                    info_dict[f"INFO_{kinfo}"] = val
        except Exception:
            pass

        row = {
            "vcf_source": sample_name,
            "chrom": chrom,
            "pos": pos1,
            "ref": ref,
            "alt": alt,
            "variant_type": var_type,
            "qual": qual,
            "filter": flt,
            "kmer": kmer,
            "gc_window": gc,
            "region_type": region,
            "gene_id": gene_id,
            "gene_biotype": gene_biotype,
            "bed_annot": bed_annot,
            "repeat_class": repeat_class,
            "flagger_state": flagger_state,
        }
        # extend with INFO entries
        row.update(info_dict)

        rows.append(row)

    df = pd.DataFrame(rows)
    return df

# ------------------------------------------------------------
# Preprocessing for CatBoost
# ------------------------------------------------------------
def prepare_for_catboost(df, label_col="TRUE_VCF"):
    """
    - Ensure label_col exists; otherwise error.
    - Fill missing values:
       - categorical: fill with 'NA'
       - numerical: fill with median
    - Convert categorical columns to strings (object) for CatBoost
    - Return X (DataFrame), y (Series), list_of_cat_feature_names
    """
    if label_col not in df.columns:
        raise ValueError(f"Label column '{label_col}' not found in DataFrame. Available columns: {df.columns.tolist()}")

    # drop columns that are clearly identifiers/should not be used as features (optional)
    # keep chrom/pos maybe as group features; but we won't use 'vcf_source' as numeric
    df = df.copy()

    y = df[label_col].astype(int)
    X = df.drop(columns=[label_col])

    # Detect categorical: consider object dtype + a few known columns
    # Also treat ref/alt/variant_type/region_type/gene_biotype/bed_annot/repeat_class/flagger_state as categorical
    forced_cat = ["ref", "alt", "variant_type", "region_type", "gene_biotype", "bed_annot", "repeat_class", "flagger_state", "vcf_source", "filter"]
    cat_cols = [c for c in forced_cat if c in X.columns]
    # add any object dtype columns
    obj_cols = X.select_dtypes(include=["object"]).columns.tolist()
    for c in obj_cols:
        if c not in cat_cols:
            cat_cols.append(c)

    # numeric columns
    num_cols = [c for c in X.columns if c not in cat_cols]

    # Fill missing
    for c in cat_cols:
        X[c] = X[c].fillna("NA").astype(str)

    for c in num_cols:
        # convert to numeric where possible
        X[c] = pd.to_numeric(X[c], errors="coerce")
        median = X[c].median(skipna=True)
        X[c] = X[c].fillna(median)

    # Ensure ordering stable
    X = X.reindex(sorted(X.columns), axis=1)

    # cat feature indices for Pool
    cat_feature_indices = [i for i, c in enumerate(X.columns) if c in cat_cols]
    cat_feature_names = [c for c in X.columns if c in cat_cols]

    return X, y, cat_feature_indices, cat_feature_names

# ------------------------------------------------------------
# Training loop (GroupKFold by chrom)
# ------------------------------------------------------------
def train_catboost(X, y, cat_feature_indices, groups, model_out, n_splits=5):
    """
    Train CatBoost with GroupKFold (groups typically chromosome) to avoid leakage.
    Saves model to model_out.
    Returns trained model and aggregated predictions on test folds.
    """
    gkf = GroupKFold(n_splits=n_splits)
    oof_preds = np.zeros(len(X))
    fold = 0
    metrics = []
    for train_idx, val_idx in gkf.split(X, y, groups):
        fold += 1
        print(f"[INFO] Fold {fold} train {len(train_idx)} val {len(val_idx)}")
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        train_pool = Pool(X_train, y_train, cat_features=cat_feature_indices)
        val_pool = Pool(X_val, y_val, cat_features=cat_feature_indices)

        # class weights - inverse frequency
        pos = y_train.sum()
        neg = len(y_train) - pos
        if pos == 0:
            class_weights = None
        else:
            w_pos = neg / pos
            class_weights = [1.0, float(w_pos)]
            print(f"[INFO] class_weights: {class_weights}")

        model = CatBoostClassifier(
            iterations=2000,
            learning_rate=0.05,
            depth=6,
            eval_metric='AUC',
            loss_function='Logloss',
            early_stopping_rounds=100,
            random_seed=42,
            verbose=100,
            class_weights=class_weights
        )

        model.fit(train_pool, eval_set=val_pool, use_best_model=True)

        preds = model.predict_proba(X_val)[:, 1]
        oof_preds[val_idx] = preds

        # metrics for this fold
        precision, recall, _ = precision_recall_curve(y_val, preds)
        pr_auc = auc(recall, precision)
        roc = roc_auc_score(y_val, preds) if len(np.unique(y_val)) > 1 else np.nan
        print(f"[INFO] Fold {fold} PR-AUC: {pr_auc:.4f}, ROC-AUC: {roc:.4f}")
        metrics.append({"fold": fold, "pr_auc": pr_auc, "roc_auc": roc})

    # retrain final model on full data
    full_pool = Pool(X, y, cat_features=cat_feature_indices)
    final_model = CatBoostClassifier(
        iterations= model.get_best_iteration() if hasattr(model, "get_best_iteration") else 1000,
        learning_rate=0.05,
        depth=6,
        eval_metric='AUC',
        loss_function='Logloss',
        random_seed=42,
        verbose=100,
        class_weights=class_weights
    )
    final_model.fit(full_pool, use_best_model=False)
    final_model.save_model(str(model_out))
    print(f"[INFO] Saved model to {model_out}")

    return final_model, oof_preds, metrics

# ------------------------------------------------------------
# Main CLI
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="VCF -> feature table -> CatBoost training")
    parser.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    parser.add_argument("--gff3", required=True, help="GFF3 annotation")
    parser.add_argument("--vcf", nargs="+", required=True, help="VCF files")
    parser.add_argument("--bed", help="Generic BED (genes/other)", default=None)
    parser.add_argument("--repeat", help="RepeatMasker BED", default=None)
    parser.add_argument("--flagger", help="Flagger BED", default=None)
    parser.add_argument("--out", required=True, help="Output TSV/CSV for features")
    parser.add_argument("--model", required=True, help="Output CatBoost model file")
    parser.add_argument("--label", default="TRUE_VCF", help="Label column name inside VCF INFO")
    args = parser.parse_args()

    # load resources
    print("[INFO] Loading reference...")
    fasta = Fasta(args.ref)

    print("[INFO] Loading GFF3 DB...")
    gff_db = load_gff_db(args.gff3)

    # load bed files with pybedtools if provided
    bed_files = {}
    def _load_bed(p):
        if p is None:
            return None
        return pybedtools.BedTool(str(p))

    bed_files["genes"] = _load_bed(args.bed)
    bed_files["repeat"] = _load_bed(args.repeat)
    bed_files["flagger"] = _load_bed(args.flagger)

    # iterate VCFs
    dfs = []
    for v in args.vcf:
        print(f"[INFO] Processing VCF {v} ...")
        df_v = extract_from_vcf(v, fasta, gff_db, bed_files=bed_files)
        dfs.append(df_v)
    if len(dfs) == 0:
        print("[ERROR] no VCFs processed.")
        sys.exit(1)
    df_all = pd.concat(dfs, ignore_index=True, sort=False)

    # If label is present in INFO_... columns, standardize name
    info_label_col = None
    # possible keys: TRUE_VCF (direct column) or INFO_TRUE_VCF (from INFO)
    if args.label in df_all.columns:
        info_label_col = args.label
    elif ("INFO_" + args.label) in df_all.columns:
        info_label_col = "INFO_" + args.label
    else:
        # attempt to find any column that looks like the label ignoring case
        for c in df_all.columns:
            if c.upper().endswith(args.label.upper()):
                info_label_col = c
                break

    if info_label_col is None:
        print(f"[ERROR] Label column '{args.label}' not found in VCF INFO fields. Columns: {df_all.columns.tolist()[:50]}")
        print("[HINT] Ensure your VCF INFO contains a TRUE_VCF flag or specify --label accordingly.")
        sys.exit(1)

    # rename label column to standard name
    df_all = df_all.rename(columns={info_label_col: "TRUE_VCF"})

    # save features raw
    df_all.to_csv(args.out, sep="\t", index=False)
    print(f"[INFO] Raw features saved to {args.out}")

    # prepare for catboost
    print("[INFO] Preparing data for CatBoost...")
    X, y, cat_idxs, cat_names = prepare_for_catboost(df_all, label_col="TRUE_VCF")
    print(f"[INFO] cat features ({len(cat_names)}): {cat_names}")

    # groups for GroupKFold - use chrom
    groups = X["chrom"].astype(str).values if "chrom" in X.columns else np.arange(len(X))

    # train
    print("[INFO] Training CatBoost with GroupKFold by chromosome...")
    model_out = Path(args.model)
    model, oof_preds, metrics = train_catboost(X, y, cat_idxs, groups, model_out, n_splits=5)

    # attach predictions and save
    df_all["pred_prob"] = oof_preds
    preds_out = str(model_out.with_suffix(".preds.csv"))
    df_all.to_csv(preds_out, sep="\t", index=False)
    print(f"[INFO] Predictions saved to {preds_out}")

if __name__ == "__main__":
    main()
