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
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
# ------------------------------------------------------------
# SVTYPE normalization
# ------------------------------------------------------------

SV_TYPES = {"DEL","INS","INV","DUP","BND"}

def normalize_svtype(x):

    if x in SV_TYPES:
        return "SV"

    return "SNV"

# ------------------------------------------------------------
# CatBoost
# ------------------------------------------------------------

from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, precision_recall_curve
from sklearn.calibration import CalibratedClassifierCV
import numpy as np


def train_catboost(df, label="TRUE_VARIANT", alpha=5.0):

    drop_cols = {
        label, "CHROM", "POS", "END",
        "REF", "ALT", "SVTYPE", "SVTYPE_NORM"
    }

    X = df.drop(columns=[c for c in drop_cols if c in df.columns])
    y = df[label]

    if y.nunique() < 2:
        raise RuntimeError("Only one class present")

    cat_cols = X.select_dtypes(include="object").columns.tolist()

    # 🔬 split: train / val / test
    X_tr, X_tmp, y_tr, y_tmp = train_test_split(
        X, y, test_size=0.8, stratify=y, random_state=42
    )

    X_val, X_te, y_val, y_te = train_test_split(
        X_tmp, y_tmp, test_size=0.5, stratify=y_tmp, random_state=42
    )

    # ансамбль с разными параметрами
    # param_grid = [
    #     dict(depth=6, lr=0.03, l2=10, bt=1, rsm=0.8),
    #     dict(depth=8, lr=0.02, l2=20, bt=2, rsm=0.7),
    #     dict(depth=4, lr=0.05, l2=5,  bt=0, rsm=1.0),
    #     dict(depth=10,lr=0.03, l2=30, bt=1, rsm=0.6),
    #     dict(depth=6, lr=0.06, l2=40, bt=3, rsm=0.9),
    # ]

    # models = []

    # for i, p in enumerate(param_grid):
    #     print(f"Training model {i+1}")

    #     model = CatBoostClassifier(
    #         iterations=600,
    #         depth=p["depth"],
    #         learning_rate=p["lr"],
    #         l2_leaf_reg=p["l2"],
    #         bagging_temperature=p["bt"],
    #         rsm=p["rsm"],
    #         loss_function="Focal:focal_alpha=0.25;focal_gamma=3",
    #         eval_metric="PRAUC",
    #         # auto_class_weights="SqrtBalanced",   # усиливаем precision
    #         random_seed=42 + i,
    #         verbose=100
        # )
    # ансамбль по разным сидов
    param_grid = [
        dict(random_seed=11, depth=2),
        dict(random_seed=22, depth=3),
        dict(random_seed=33, depth=4),
        dict(random_seed=44, depth=4),
        dict(random_seed=55, depth=6),
        dict(random_seed=23, depth=3),
        dict(random_seed=42, depth=8),
        dict(random_seed=43, depth=4),
        dict(random_seed=58, depth=2),
        dict(random_seed=102, depth=6),
        dict(random_seed=105, depth=2),
        dict(random_seed=360, depth=3),
        dict(random_seed=1054, depth=4),
        dict(random_seed=78, depth=4),
        dict(random_seed=5557, depth=6),
        dict(random_seed=455, depth=3),
        dict(random_seed=1258, depth=8),
        dict(random_seed=657, depth=4),
        dict(random_seed=1257, depth=2),
        dict(random_seed=4568, depth=6),
    ]

    models = []
    calibrated_models = []

    for i, params in enumerate(param_grid):
        print(f"\nTraining model {i+1}")

        model = CatBoostClassifier(
            iterations=300,
            learning_rate=0.05,
            loss_function="Focal:focal_alpha=0.25;focal_gamma=3",
            l2_leaf_reg=50,
            subsample=0.7,
            rsm=0.8,
            random_strength=2,
            bagging_temperature=2,
            eval_metric="PRAUC",
            verbose=100,
            posterior_sampling=True,
            **params
        )

        model.fit(
            X_tr,
            y_tr,
            cat_features=cat_cols,
            eval_set=(X_val, y_val)
        )

        models.append(model)

    # ensemble prediction
    probs = np.mean(
        [m.predict_proba(X_te)[:, 1] for m in models],
        axis=0
    )

    auc = roc_auc_score(y_te, probs)
    print("\nEnsemble AUC =", round(auc, 4))

    fi = model.get_feature_importance(prettified=True)
    print(fi.head(5))

    val_probs = np.mean(
        [m.predict_proba(X_val)[:, 1] for m in models],
        axis=0
    )
    
    precision, recall, thresholds = precision_recall_curve(y_val, val_probs)

    mask = precision[:-1] >= 0.95

    best_idx = np.argmax(recall[:-1][mask])
    best_threshold = thresholds[mask][best_idx]
    best_recall = recall[:-1][mask][best_idx]

    print(f"Best threshold: {best_threshold:.4f}")
    print(f"Recall @ precision≥0.95: {best_recall:.4f}")

    return {
        "models": models,
        "threshold": best_threshold,
        "cat_cols": cat_cols
    }

def prepare_X(X, cat_cols):
    X = X.copy()

    for c in cat_cols:
        if c in X.columns:
            X[c] = X[c].astype(str)

    return X

def predict_block(model_dict, X, label="TRUE_VARIANT"):

    drop_cols = {
        label, "CHROM", "POS", "END",
        "REF", "ALT", "SVTYPE", "SVTYPE_NORM"
    }

    X = X.drop(columns=[c for c in drop_cols if c in X.columns])

    probs = np.mean(
        [m.predict_proba(X)[:, 1] for m in model_dict["models"]],
        axis=0
    )

    return probs


# ------------------------------------------------------------
# VCF writing
# ------------------------------------------------------------

def write_vcf(df,vcfs,out_vcf,threshold):

    df=df[df.ML_PROB>=threshold]

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

def fix_categories(df):
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)
    return df

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main(args):

    print("[1] Loading TSV")
    df = pd.read_csv(args.tsv, sep="\t")
    print(f'{args.prefix}')


    print("[2] Splitting SNV / SV")
    df_snv = df[df.SVTYPE_NORM == "SNV"].copy()
    df_snv = fix_categories(df_snv)
    # df_sv  = df[df.SVTYPE_NORM == "SV"].copy()

    print("SNV:", len(df_snv))
    # print("SV :", len(df_sv))

    # --------------------------------------------------------
    # Models to try
    # --------------------------------------------------------

    # models_to_try = [
    #     {
    #         "name": "focal",
    #         "params": {
    #             "iterations": 300,
    #             "depth": 10,
    #             "learning_rate": 0.03,
    #             "loss_function":"'Focal:focal_alpha=0.75;focal_gamma=3'",
    #             "auto_class_weights": "Balanced",
    #             "l2_leaf_reg": 8,
    #             "eval_metric": "AUC",
    #             "verbose": 100
    #         }
    #     },
    #     {
    #         "name": "light",
    #         "params": {
    #             "iterations": 300,
    #             "depth": 8,
    #             "learning_rate": 0.1,
    #             "loss_function": "Logloss",
    #             "auto_class_weights": "Balanced",
    #             "l2_leaf_reg": 3,
    #             "eval_metric": "AUC",
    #             "verbose": 100
    #         }
    #     }
    # ]

    # --------------------------------------------------------
    # Run all models
    # --------------------------------------------------------

    # for model_cfg in models_to_try:

        # df_snv = df_snv.drop(columns=["ML_PROB"])
        # name = model_cfg["name"]
        # params = model_cfg["params"]

        # print(f"\n=== Training model: {name} ===")

        # # SNV model
        # model_snv = train_catboost(df_snv, params)
        # # SV model
        # # model_sv  = train_catboost(df_sv,  params)

        # print(f"[7] prediction for {name}")

        # df_snv[f"ML_PROB"] = 0.0
        # df_snv.loc[df_snv.index, "ML_PROB"] = predict_block(df_snv, model_snv)
        # # df.loc[df_sv.index,  f"ML_PROB_{name}"] = predict_block(df_sv,  model_sv)

        # # Save TSV
        # out_tsv = f"{name}.variants.tsv"
        # df_snv.to_csv(out_tsv, sep="\t", index=False)
        # print(f"Saved {out_tsv}")

        # # Save VCFs
        # print(f"[8] VCF writing for {name}")

        # df_snv_sub = df_snv[df_snv.SVTYPE_NORM == "SNV"]
        # # df_sv_sub  = df[df.SVTYPE_NORM == "SV"]

        # write_vcf(df_snv_sub, args.vcfs, f"{name}.snv.vcf", args.threshold)
        # # write_vcf(df_sv_sub,  args.vcfs, f"{name}.sv.vcf",  args.threshold)
        # # write_vcf(df,         args.vcfs, f"{name}.all.vcf", args.threshold)

        # print(f"Model {name} done.")
    
    df_snv = df_snv.drop(columns=["ML_PROB"])
    df_snv[f"ML_PROB"] = 0.0
    print(f"\n=== Training model")

    # SNV model
    # model_snv = train_catboost(df_snv, params)
    model_snv = train_catboost(df_snv)
    # SV model
    # model_sv  = train_catboost(df_sv,  params)

    print(f"[7] prediction for model")

    df_snv.loc[df_snv.index, "ML_PROB"] = predict_block(model_snv, df_snv)
    # df.loc[df_sv.index,  f"ML_PROB_{name}"] = predict_block(df_sv,  model_sv)

    # Save TSV
    out_tsv = f"{args.prefix}.variants.tsv"
    df_snv.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved {out_tsv}")

    # Save VCFs
    print(f"[8] VCF writing for ensemble")

    df_snv_sub = df_snv[df_snv.SVTYPE_NORM == "SNV"]
    # df_sv_sub  = df[df.SVTYPE_NORM == "SV"]

    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}'.snv.vcf", model_snv["threshold"])
    # write_vcf(df_sv_sub,  args.vcfs, f"{name}.sv.vcf",  args.threshold)
    # write_vcf(df,         args.vcfs, f"{name}.all.vcf", args.threshold)

    print(f"Model ensemble done.")

    df = df_snv
    plt.hist(df.loc[df.TRUE_VARIANT == 0, "ML_PROB"], bins=50, alpha=0.5, label="0")
    plt.hist(df.loc[df.TRUE_VARIANT == 1, "ML_PROB"], bins=50, alpha=0.5, label="1")
    plt.legend()
    plt.yscale("log")

    plt.savefig("hist_seed.png", dpi=200)

    plt.figure()
    y_true = df.TRUE_VARIANT.values
    y_score = df.ML_PROB.values

    prec, rec, thr = precision_recall_curve(y_true, y_score)

    plt.plot(rec, prec)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.savefig("pr_curve_seed.png", dpi=200)

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True, help="Input TSV from main pipeline")
    p.add_argument("--vcfs",nargs="+",required=True)
    p.add_argument("--prefix",required=True)
    p.add_argument("--threshold",type=float,default=0.9)
    main(p.parse_args())

