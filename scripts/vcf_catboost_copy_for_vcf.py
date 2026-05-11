#!/usr/bin/env python3

import argparse
import pysam
from pathlib import Path
import pandas as pd
import numpy as np
from intervaltree import IntervalTree
from cyvcf2 import VCF
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.calibration import calibration_curve
from sklearn.metrics import roc_curve
import random
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from functools import reduce


# ------------------------------------------------------------
# utilities
# ------------------------------------------------------------
def fix_categories(df):
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)
    return df

def save_hyperparam_results(results, path="hyperparam_results.tsv"):
    """
    Сохраняет результаты подбора гиперпараметров в TSV.
    results — список словарей вида:
        {
            "params": {...},
            "recall": float,
            "threshold": float,
            "model": CatBoostClassifier
        }
    """
    rows = []

    for r in results:
        row = {
            "recall": r["recall"],
            "threshold": r["threshold"]
        }

        # разворачиваем словарь params в отдельные колонки
        for k, v in r["params"].items():
            row[k] = v

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    print(f"Saved hyperparameter search results to {path}")

def safe_float(x):
    try:
        return float(x)
    except:
        return 0.0


def normalize_info_value(key, value):

    if value is None:
        return 0.0

    if key.upper() == "SVLEN":
        if isinstance(value, (list, tuple)):
            return abs(safe_float(value[0]))
        return abs(safe_float(value))

    if isinstance(value, (list, tuple)):
        try:
            return float(np.mean(value))
        except:
            return 0.0

    return safe_float(value)


# ------------------------------------------------------------
# SVTYPE normalization
# ------------------------------------------------------------

SV_TYPES = {"DEL","INS","INV","DUP","BND"}

def normalize_svtype(x):

    if x in SV_TYPES:
        return "SV"

    return "SNV"


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

        if fmt == "GT" and isinstance(val, (list,np.ndarray)):
            out[f"{prefix}_GT"] = int(val[0] != val[1])

        elif isinstance(val,(int,float,np.integer,np.floating)):
            out[key] = float(val)

        elif isinstance(val,(list,np.ndarray)):
            nums = [x for x in val if isinstance(x,(int,float))]
            out[key] = float(np.mean(nums)) if nums else 0.0

        else:
            out[key] = str(val)

    return out


# ------------------------------------------------------------
# VCF loading
# ------------------------------------------------------------
def load_telomers(bed_path):
    telomeres = {}

    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)

            if chrom not in telomeres:
                telomeres[chrom] = []

            telomeres[chrom].append((start, end))

    return telomeres

def is_excluded(chrom, pos, end, exclude_dict):
    for start, stop in exclude_dict[chrom]:
        if not (end < start or pos > stop):
            return True

    return False

def drop_useless_columns(df):
    cols_to_drop = []

    for col in df.columns:

        if col == "SVTYPE" or  col == "SVLEN" or col.endswith("_present"):
            continue
        uniq = set(df[col].astype(str).unique())

        # если только одно уникальное значение → мусор
        if len(uniq) == 1:
            cols_to_drop.append(col)

    return df.drop(columns=cols_to_drop), cols_to_drop


def load_vcf_as_df(path, exclude_trees=None):

    prefix = Path(path).stem
    rows = []

    for v in VCF(path):

        info = dict(v.INFO)
        
        if exclude_trees is not None:
            if is_excluded(v.CHROM, v.POS, v.END, exclude_trees):
                continue  

        
        row = {
            "CHROM": v.CHROM,
            "POS": v.POS,
            "END": max(v.end, v.POS),
            "REF": v.REF,
            "ALT": str(v.ALT[0]) if v.ALT else None,
            "SVTYPE": info.get("SVTYPE","SNV"),
            "SVLEN": normalize_info_value("SVLEN", info.get("SVLEN")),
            f"{prefix}_present":1
        }

        for k,val in info.items():
            row[f"{prefix}_INFO_{k}"] = normalize_info_value(k,val)

        row.update(parse_format_fields(v,prefix))

        rows.append(row)

    df = pd.DataFrame(rows)

    df, dropped = drop_useless_columns(df)

    if dropped:
        print(f"Удалены бесполезные столбцы: {dropped}")

    return df

def outer_merge_vcfs(paths, exclude_trees=None):

    dfs = [load_vcf_as_df(p, exclude_trees) for p in paths]

    df = dfs[0]

    for d in dfs[1:]:

        df = df.merge(
            d,
            on=["CHROM", "POS"], 
            how="outer"
        )
        for col in ["REF", "ALT", "SVTYPE", "SVLEN", "END"]:
            col_x = f"{col}_x"
            col_y = f"{col}_y"

            if col_x in df.columns and col_y in df.columns:
                df[col] = df[col_x].combine_first(df[col_y])
                df.drop(columns=[col_x, col_y], inplace=True)

    return df.fillna(0)



# ------------------------------------------------------------
# MERFIN
# ------------------------------------------------------------

def load_merfin_support(vcf):

    return {(v.CHROM,v.POS) for v in VCF(vcf)}


def annotate_merfin_support(df, merfin_vcf):

    support = load_merfin_support(merfin_vcf)

    df["MERFIN_SUPPORTED"] = [
    (c,p) in support for c,p in zip(df.CHROM,df.POS)
]

    df["MERFIN_SUPPORTED"] = df["MERFIN_SUPPORTED"].astype(int)

    return df


def load_merfin_scores(bed):

    trees = {}

    with open(bed) as f:

        for line in f:

            if line.startswith("#"):
                continue

            chrom,start,end,score = line.split()[:4]

            trees.setdefault(chrom,IntervalTree()).addi(
                int(start),
                int(end),
                float(score)
            )

    return trees


def annotate_merfin_scores(df,bed):

    trees = load_merfin_scores(bed)

    scores = np.zeros(len(df))

    for i,(chrom,pos,end) in enumerate(
        zip(df.CHROM,df.POS,df.END)
    ):

        hits = trees.get(chrom,IntervalTree()).overlap(pos,end)

        if hits:
            scores[i] = next(iter(hits)).data

    df["MERFIN_SCORE"] = scores

    return df


# ------------------------------------------------------------
# annotations
# ------------------------------------------------------------

def load_bed(path,value_col=None):

    trees={}

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            parts=line.rstrip().split("\t")

            if len(parts)<3:
                continue

            chrom,start,end = parts[:3]

            if int(start) >= int(end):
                continue

            val = parts[value_col] if value_col else 1

            trees.setdefault(chrom,IntervalTree()).addi(
                int(start),
                int(end),
                val
            )

    return trees


def load_gff3(path,attr_key=None):

    trees={}

    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue

            parts=line.rstrip().split("\t")

            if len(parts)<9:
                continue

            chrom=parts[0]
            start=int(parts[3])
            end=int(parts[4])

            if int(start) >= int(end):
                continue


            value=parts[2]

            if attr_key:

                attrs=dict(
                    x.split("=",1)
                    for x in parts[8].split(";")
                    if "=" in x
                )

                raw=attrs.get(attr_key)

                if raw:
                    value=raw.split()[0]

            trees.setdefault(chrom,IntervalTree()).addi(
                start,end,value
            )

    return trees


def annotate_interval(df,trees,col):

    present=np.zeros(len(df))
    value=["NA"]*len(df)

    for i,(chrom,pos,end) in enumerate(
        zip(df.CHROM,df.POS,df.END)
    ):

        hits=trees.get(chrom,IntervalTree()).overlap(pos,end)

        if hits:

            h=next(iter(hits))

            present[i]=1
            value[i]=h.data

    df[f"{col}_present"]=present
    df[col]=value


# ------------------------------------------------------------
# truth
# ------------------------------------------------------------

def load_truth(vcf):
    truth = {}
    opener = gzip.open if vcf.endswith(".gz") else open

    with opener(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, *_ = line.split("\t")
            truth.setdefault(chrom, []).append(int(pos))

    for k in truth:
        truth[k] = np.array(truth[k], dtype=np.int32)

    return truth


def build_truth_labels(df, truth):
    labels = np.zeros(len(df), dtype=int)

    for chrom, idx in df.groupby("CHROM").groups.items():
        if chrom not in truth:
            continue

        pos = np.sort(df.loc[idx, "POS"].values)
        t = truth[chrom]

        i = j = 0
        hits = np.zeros(len(pos), dtype=bool)

        while i < len(pos) and j < len(t):
            if int(pos[i]) == int(t[j]):
                hits[i] = True
                i += 1
            elif t[j] < pos[i] - 1:
                j += 1
            else:
                i += 1

        # вернуть в исходный порядок
        labels[idx] = hits[np.argsort(np.argsort(df.loc[idx, "POS"].values))]

    return labels


# ------------------------------------------------------------
# CatBoost
# ------------------------------------------------------------

import pandas as pd

def downsample_negatives(df, target_col="TRUE_VARIANT", neg_frac=None, neg_to_pos_ratio=None, random_state=42):
    """
    Downsample negative class.
    
    Parameters:
        df: DataFrame
        target_col: name of target column
        neg_frac: fraction of negatives to keep (e.g., 0.1)
        neg_to_pos_ratio: keep N negatives per 1 positive (e.g., 5)
    """
    pos = df[df[target_col] == 1]
    neg = df[df[target_col] == 0]

    if neg_frac is not None:
        neg_down = neg.sample(frac=neg_frac, random_state=random_state)
    elif neg_to_pos_ratio is not None:
        n = min(len(neg), len(pos) * neg_to_pos_ratio)
        neg_down = neg.sample(n=n, random_state=random_state)
    else:
        raise ValueError("Specify either neg_frac or neg_to_pos_ratio")

    df_balanced = pd.concat([pos, neg_down]).sample(frac=1, random_state=random_state)
    return df_balanced

def hard_negative_mining(model, df, features, target_col="TRUE_VARIANT",
                         threshold=None, top_k=None, ratio=None):
    """
    Hard Negative Mining for CatBoost (or any classifier with predict_proba).
    
    Parameters:
        model: trained model
        df: full dataset
        features: list of feature columns
        threshold: keep negatives with pred > threshold
        top_k: keep top K hardest negatives
        ratio: keep ratio * number_of_positives negatives
    """
    pos = df[df[target_col] == 1]
    neg = df[df[target_col] == 0].copy()

    # predict probabilities for negatives
    neg["pred"] = model.predict_proba(neg[features])[:, 1]

    if threshold is not None:
        hard_neg = neg[neg["pred"] > threshold]

    elif top_k is not None:
        hard_neg = neg.sort_values("pred", ascending=False).head(top_k)

    elif ratio is not None:
        k = int(len(pos) * ratio)
        hard_neg = neg.sort_values("pred", ascending=False).head(k)

    else:
        raise ValueError("Specify threshold, top_k, or ratio")

    # combine positives + hard negatives
    df_hnm = pd.concat([pos, hard_neg]).sample(frac=1, random_state=42)
    return df_hnm

# ------------------------------------------------------------
# VCF writing
# ------------------------------------------------------------

def write_vcf(df,vcfs,out_vcf,threshold,param):

    df=df[df[param] >=threshold]

    if df.empty:
        return

    df=df.sort_values("ML_PROB",ascending=False)
    df=df.drop_duplicates(["CHROM","POS"])

    best={(r.CHROM,r.POS):r for r in df.itertuples()}

    template=vcfs[0]

    header=pysam.VariantFile(template).header.copy()

    if "ML_PROB" not in header.info:

        header.info.add(
            "ML_PROB",
            1,
            "Float",
            "CatBoost probability"
        )

    out=pysam.VariantFile(out_vcf,"w",header=header)

    for vcf in vcfs:

        vf=pysam.VariantFile(vcf)

        for v in vf.fetch():

            key=(v.chrom,v.pos)

            if key not in best:
                continue

            r=best[key]

            nv=out.new_record(
                contig=v.chrom,
                start=v.start,
                stop=v.stop,
                alleles=v.alleles
            )

            nv.info["ML_PROB"]=float(r.ML_PROB)

            out.write(nv)

    out.close()


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main(args):

    print("[1] loading VCFs")

    df = pd.read_csv('/mnt/data/vanichkinao/MetaPolisher/auto_pipeline/second_try/try_wh_merqury_and_logloss.variants.tsv', sep="\t")
    df = df.fillna(0)

    df["ML_PROB"] = df ["ML_PROB_one"]
    df["ML_PRED"] = df ["ML_PRED_one"]
    
    print(f"[8] VCF writing for ensemble")

    df_snv_sub = df[df.SVTYPE_NORM == "SNV"]



    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}.one_model.vcf", top_models[1]["threshold"], "ML_PROB")

    top_models[1]["model"].save_model("model_{args.prefix}.cbm")

    print(f"Model ensemble done.")

    plt.hist(df.loc[df.TRUE_VARIANT == 0, "ML_PROB"], bins=50, alpha=0.5, label="0")
    plt.hist(df.loc[df.TRUE_VARIANT == 1, "ML_PROB"], bins=50, alpha=0.5, label="1")
    plt.legend()
    plt.yscale("log")

    plt.savefig(f"hist_{args.prefix}.png", dpi=200)

    plt.figure()
    y_true = df.TRUE_VARIANT.values
    y_score = df.ML_PROB.values

    prec, rec, thr = precision_recall_curve(y_true, y_score)

    plt.plot(rec, prec)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.savefig(f"pr_curve_{args.prefix}.png", dpi=200)

    plt.figure()
    plt.hist(df.loc[df.TRUE_VARIANT == 0, "ML_PROB"], bins=50, alpha=0.5, label="0")
    plt.hist(df.loc[df.TRUE_VARIANT == 1, "ML_PROB"], bins=50, alpha=0.5, label="1")
    plt.legend()
    plt.yscale("log")

    plt.savefig(f"hist_wo_log_{args.prefix}.png", dpi=200)

    # y_val — истинные метки
    # model — обученный CatBoostClassifier
    # X_val — валидационные признаки

    y_pred_proba = top_models[1]["model"].predict_proba(X_val)[:, 1]

    # -----------------------------
    # 1. ROC Curve
    # -----------------------------
    fpr, tpr, _ = roc_curve(y_val, y_pred_proba)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"roc_curve_{args.prefix}.png", dpi=200)

    # -----------------------------
    # 2. Precision–Recall Curve
    # -----------------------------
    precision, recall, _ = precision_recall_curve(y_val, y_pred_proba)

    plt.figure(figsize=(6, 5))
    plt.plot(recall, precision)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision–Recall Curve")
    plt.grid(True)
    plt.savefig(f"precision_recall_curve_{args.prefix}.png", dpi=200)

    # -----------------------------
    # 3. Calibration Curve
    # -----------------------------
    prob_true, prob_pred = calibration_curve(y_val, y_pred_proba, n_bins=10)

    plt.figure(figsize=(6, 5))
    plt.plot(prob_pred, prob_true, marker="o")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("Predicted probability")
    plt.ylabel("True probability")
    plt.title("Calibration Curve")
    plt.grid(True)
    plt.savefig(f"calibration_curve_{args.prefix}.png", dpi=200)


    # Confusion Matrix
    cm = confusion_matrix(df["TRUE_VARIANT"], df["ML_PRED_one)"])
    labels = ["Negative", "Positive"]

    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=labels, yticklabels=labels)
    plt.xlabel("Predicted label")
    plt.ylabel("True label")
    plt.title(f"Confusion Matrix (threshold = {threshold})")
    plt.savefig(f"matrix_{args.prefix}.png", dpi=200)


if __name__=="__main__":

    p=argparse.ArgumentParser()

    p.add_argument("--vcfs",nargs="+",required=True)
    p.add_argument("--truth_vcf",required=True)
    p.add_argument("--merfin_pass_vcf")
    p.add_argument("--merfin_scores_bed")
    p.add_argument("--merqury")
    p.add_argument("--quast")
    p.add_argument("--repeat_gff")
    p.add_argument("--liftoff_gff")
    p.add_argument("--low_complex")
    p.add_argument("--flagger",required=True)
    p.add_argument("--out_vcf")
    p.add_argument("--out_table",default="variant_features.tsv")
    p.add_argument("--exclude_telomers", help="BED file with telomere regions to exclude")
    p.add_argument("--prefix",required=True)


    main(p.parse_args())