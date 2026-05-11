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
import random
# ------------------------------------------------------------
# SVTYPE normalization
# ------------------------------------------------------------

SV_TYPES = {"DEL","INS","INV","DUP","BND"}
import pandas as pd

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
            loss_function="Logloss",
            auto_class_weights="SqrtBalanced",
            l2_leaf_reg=30,
            eval_metric="AUC",
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

    return {
        "models": models,
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
    df = df.fillna(str(0.0))


    # Какие хромосомы реально есть
    chroms_present = [c for c in df["CHROM"].unique()]

    # Разбиение
    train_chroms = chroms_present[:10]      # chr1–chr16
    val_chroms   = chroms_present[10:17]    # chr17–chr19
    test_chroms  = chroms_present[17:]      # chr20–chr22, X, Y, M

    train_df = df[df["CHROM"].isin(train_chroms)]
    val_df   = df[df["CHROM"].isin(val_chroms)]
    test_df  = df[df["CHROM"].isin(test_chroms)]

    print(len(train_df), len(val_df), len(test_df))

    print("[2] Splitting SNV / SV")
    
    train_df = train_df.drop(columns=["ML_PROB"])
    train_df[f"ML_PROB"] = 0.0
    test_df= test_df.drop(columns=["ML_PROB"])
    test_df[f"ML_PROB"] = 0.0
    val_df= val_df.drop(columns=["ML_PROB"])
    val_df[f"ML_PROB"] = 0.0

    features = [c for c in df.columns if c not in ["TRUE_VARIANT", "CHROM", "POS", "END", "REF", "ALT", "SVTYPE", "SVTYPE_NORM"]]

    print(f"\n=== Training model")

    X_train = train_df[features]
    y_train = train_df["TRUE_VARIANT"]
    X_val   = val_df[features]
    y_val   = val_df["TRUE_VARIANT"]
    print(len(X_train), len(y_train), len(X_val), len(y_val))
    # -----------------------------
    # 1. Гиперпараметры для перебора
    # -----------------------------
    param_grid = []

    loss_functions = ["Logloss"]
    l2_values = [3, 10, 30]
    auto_class = ["Balanced", "SqrtBalanced"]
    learning_rates = [0.03, 0.05, 0.1]
    depths = [3, 4, 6, 8]

    for loss in loss_functions:
        for l2 in l2_values:
            for ac in auto_class:
                for lr in learning_rates:
                    for d in depths:
                            param_grid.append(dict(
                                loss_function=loss,
                                l2_leaf_reg=l2,
                                auto_class_weights=ac,
                                learning_rate=lr,
                                depth=d,
                                random_seed=random.randint(0, 10**9),
                                thread_count=8
                            ))

    print("Всего комбинаций:", len(param_grid))

    # -----------------------------
    # 2. Обучение и валидация
    # -----------------------------
    results = []
    models = []

    cat_cols_2 = ["chm13_v0.9_FORMAT_C", "flagger", "chm13_v0.9_FORMAT_GT", "deepvariant_merged_hifi_illumina_FORMAT_GT"]
    cat_cols = X_train.select_dtypes(include="object").columns.tolist() + cat_cols_2
    
    print(cat_cols)

    for i, params in enumerate(param_grid):
        print(f"\n[{i+1}/{len(param_grid)}] Обучение модели с параметрами: {params}")

        model = CatBoostClassifier(
            iterations=300,
            eval_metric="PRAUC",
            verbose=False,
            posterior_sampling=True,
            **params,
            used_ram_limit='100gb'
        )

        model.fit(X_train, y_train, eval_set=(X_val, y_val), cat_features=cat_cols, use_best_model=True)

        # вероятности на валидации
        val_prob = model.predict_proba(X_val)[:, 1]

        # подбор threshold под precision >= 0.9
        precision, recall, thresholds = precision_recall_curve(y_val, val_prob)

        target_precision = 0.90
        best_recall = 0
        best_threshold = 0

        for p, r, t in zip(precision, recall, np.append(thresholds, 1.0)):
            if p >= target_precision and r > best_recall:
                best_recall = r
                best_threshold = t

        results.append({
            "params": params,
            "model": model,
            "recall": best_recall,
            "threshold": best_threshold
        })
    
    save_hyperparam_results(results, "ensemble_hyperparams.tsv")

    N = 10
    top_models = sorted(results, key=lambda x: x["recall"], reverse=True)[:N]

    val_probs_list = [m["model"].predict_proba(X_val)[:, 1] for m in top_models]
    val_probs_ens = np.mean(val_probs_list, axis=0)

    precision, recall, thresholds = precision_recall_curve(y_val, val_probs_ens)

    best_recall = 0
    best_threshold = 0

    for p, r, t in zip(precision, recall, np.append(thresholds, 1.0)):
        if p >= 0.90 and r > best_recall:
            best_recall = r
            best_threshold = t

    print("\nАнсамбль:")
    print("Recall:", best_recall)
    print("Threshold:", best_threshold)


    print(f"[7] prediction for ensemble")

    # вероятности всех моделей на test
    test_probs_list = [
        m.predict_proba(test_df[features])[:, 1]
        for m in top_models["model"]
    ]

    # усреднение
    test_prob = np.mean(test_probs_list, axis=0)

    # бинарные предсказания
    test_pred = (test_prob >= best_threshold).astype(int)

    from sklearn.metrics import classification_report
    print(classification_report(test_df["TRUE_VARIANT"], test_pred))


    # вероятности ансамбля по всему геному
    all_probs_list = [
        m.predict_proba(df[features])[:, 1]
        for m in top_models
    ]

    df["ML_PROB"] = np.mean(all_probs_list, axis=0)
    df["ML_PRED"] = (df["ML_PROB"] >= best_threshold).astype(int)

   
    out_tsv = f"{args.prefix}.variants.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved {out_tsv}")

    # Save VCFs
    print(f"[8] VCF writing for ensemble")

    df_snv_sub = df[df_snv.SVTYPE_NORM == "SNV"]

    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}'.snv.vcf", best_threshold)

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

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True, help="Input TSV from main pipeline")
    p.add_argument("--vcfs",nargs="+",required=True)
    p.add_argument("--prefix",required=True)
    p.add_argument("--threshold",type=float,default=0.9)
    main(p.parse_args())
