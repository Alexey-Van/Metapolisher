#!/usr/bin/env python3

import argparse
import pysam
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
from intervaltree import IntervalTree
from cyvcf2 import VCF
from catboost import CatBoostClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.calibration import calibration_curve
from sklearn.metrics import confusion_matrix
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

    if rows:
        df = pd.DataFrame(rows)
    else:
        df = pd.DataFrame(columns=[
            "recall", "threshold"
        ])
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
    if chrom != 'chrM':
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
            if is_excluded(v.CHROM, v.POS, v.end, exclude_trees):
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

    if df.empty:
        return pd.DataFrame({
            "CHROM": pd.Series(dtype="object"),
            "POS": pd.Series(dtype="int64"),
            "END": pd.Series(dtype="int64"),
            "REF": pd.Series(dtype="object"),
            "ALT": pd.Series(dtype="object"),
            "SVTYPE": pd.Series(dtype="object"),
            "SVLEN": pd.Series(dtype="float64"),
            f"{prefix}_present": pd.Series(dtype="int64"),
        })

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

            val = parts[value_col] if value_col is not None and value_col < len(parts) else 1

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

        hits=trees.get(chrom,IntervalTree()).overlap(pos-1,end)

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

def train_catboost(df,label="TRUE_VARIANT"):

    drop_cols={
        label,"CHROM","POS","END",
        "REF","ALT","SVTYPE","SVTYPE_NORM"
    }

    X=df.drop(columns=[c for c in drop_cols if c in df.columns])
    y=df[label]

    if y.nunique()<2:
        raise RuntimeError("Only one class present")

    cat_cols = X.select_dtypes(include=["object", "string"]).columns.tolist()

    X_tr,X_te,y_tr,y_te=train_test_split(
        X,y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    model=CatBoostClassifier(
        iterations=300,
        depth=6,
        learning_rate=0.05,
        loss_function="Logloss",
        auto_class_weights="Balanced",
        l2_leaf_reg=5,
        eval_metric="AUC",
        verbose=100
    )

    model.fit(
        X_tr,
        y_tr,
        cat_features=cat_cols,
        eval_set=(X_te,y_te)
    )

    auc=roc_auc_score(y_te,model.predict_proba(X_te)[:,1])

    print("AUC =",round(auc,4))

    fi = model.get_feature_importance(prettified=True)
    print(fi.head(10))
    
    return model


def predict_block(df,model):

    drop_cols={
        "TRUE_VARIANT","CHROM","POS","END",
        "REF","ALT","SVTYPE","SVTYPE_NORM"
    }
    fi = model.get_feature_importance(prettified=True)
    print(fi.head(5))

    X=df.drop(columns=[c for c in drop_cols if c in df.columns])

    return model.predict_proba(X)[:,1]


# ------------------------------------------------------------
# VCF writing
# ------------------------------------------------------------

def write_vcf(df,vcfs,out_vcf,threshold,param):

    df=df[df[param] >=threshold]

    df=df.sort_values(param,ascending=False)
    df=df.drop_duplicates(["CHROM","POS"])

    best={(r.CHROM,r.POS):r for r in df.itertuples()}

    template=vcfs[0]
    with pysam.VariantFile(template) as template_vcf:
        header=template_vcf.header.copy()

    if param not in header.info:

        header.info.add(
            param,
            1,
            "Float",
            "CatBoost probability"
        )

    out=pysam.VariantFile(out_vcf,"w",header=header)

    emitted = set()
    for vcf in vcfs:
        with pysam.VariantFile(vcf) as vf:
            for v in vf:
                key=(v.chrom,v.pos)

                if key not in best or key in emitted:
                    continue

                r=best[key]
                emitted.add(key)

                nv=out.new_record(
                    contig=v.chrom,
                    start=v.start,
                    stop=v.stop,
                    alleles=v.alleles
                )

                nv.info[param]=float(r.ML_PROB)

                out.write(nv)

    out.close()


def predict_with_catboost(args):

    if not args.ml_model:
        raise ValueError("predict mode requires --ml_model pointing to a CatBoost .cbm file")

    if not Path(args.ml_model).is_file():
        raise FileNotFoundError(f"CatBoost model was not found: {args.ml_model}")
    if Path(args.ml_model).suffix.lower() != ".cbm":
        raise ValueError("--ml_model must point to a CatBoost .cbm file")

    print("[1] loading VCFs")
    exclude_trees = load_telomers(args.exclude_telomers) if args.exclude_telomers else None
    df = outer_merge_vcfs(args.vcfs, exclude_trees)

    print("[2] annotations")
    if args.flagger:
        annotate_interval(df, load_bed(args.flagger, 3), "flagger")
    if args.merqury:
        annotate_interval(df, load_bed(args.merqury), "merqury")

    model = CatBoostClassifier()
    model.load_model(args.ml_model)

    excluded = {
        "TRUE_VARIANT", "CHROM", "POS", "END", "REF", "ALT",
        "SVTYPE", "SVTYPE_NORM"
    }
    available = [col for col in df.columns if col not in excluded]
    expected = list(model.feature_names_ or available)
    X = df.reindex(columns=expected)

    categorical_indices = set(model.get_cat_feature_indices())
    for idx, col in enumerate(expected):
        if idx in categorical_indices:
            X[col] = X[col].fillna("0").astype(str)
        else:
            X[col] = pd.to_numeric(X[col], errors="coerce").fillna(0)

    print(f"[3] prediction with {len(expected)} model features")
    if df.empty:
        df["ML_PROB"] = pd.Series(index=df.index, dtype=float)
    else:
        df["ML_PROB"] = model.predict_proba(X)[:, 1]

    table_path = args.out_table or f"{args.prefix}.ml.variants.tsv"
    vcf_path = args.out_vcf or f"{args.prefix}.ml.vcf"
    df.to_csv(table_path, sep="\t", index=False)
    print(f"Saved {table_path}")

    # The current model's training path emits filtered SNVs. SV calls remain
    # available in the feature table and can be enabled in a future model.
    df_snv = df[df["SVTYPE"].apply(normalize_svtype) == "SNV"]
    write_vcf(df_snv, args.vcfs, vcf_path, args.threshold, "ML_PROB")
    print(f"Saved {vcf_path} (probability threshold: {args.threshold})")


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main(args):

    if args.mode == "predict":
        return predict_with_catboost(args)

    if not args.truth_vcf:
        raise ValueError("train mode requires --truth_vcf")
    if not args.flagger or not args.merqury:
        raise ValueError("train mode requires --flagger and --merqury")

    print("[1] loading VCFs")

    exclude_trees = None
    if args.exclude_telomers:
        exclude_trees = load_telomers(args.exclude_telomers)


    df=outer_merge_vcfs(args.vcfs, exclude_trees)
    if df.empty:
        raise RuntimeError("No candidate variants were found in the input VCFs")

    # print("[2] merfin")

    # df=annotate_merfin_support(df,args.merfin_pass_vcf)

    # if args.merfin_scores_bed:
    #     df=annotate_merfin_scores(df,args.merfin_scores_bed)

    print("[3] annotations")

    # repeat=load_gff3(args.repeat_gff,"Target")
    # liftoff=load_gff3(args.liftoff_gff)
    # low_complex=load_bed(args.low_complex)
    flagger=load_bed(args.flagger,3)
    merqury=load_bed(args.merqury)

    # annotate_interval(df,repeat,"repeat")
    # annotate_interval(df,liftoff,"liftoff")
    # annotate_interval(df,low_complex,"low_complex")
    annotate_interval(df,flagger,"flagger")
    annotate_interval(df,merqury,"merqury")

    print("[4] truth")

    truth=load_truth(args.truth_vcf)

    df["TRUE_VARIANT"]=build_truth_labels(df,truth)
    print("Найдено, вариантов:", sum(df["TRUE_VARIANT"]))
    class_counts = df["TRUE_VARIANT"].value_counts()
    if len(class_counts) < 2 or class_counts.min() < 3:
        raise RuntimeError(
            "Training needs at least three true and three false candidate variants "
            "to create stratified train/validation/test splits"
        )

    print("[5] SVTYPE")

    df["SVTYPE_NORM"]=df["SVTYPE"].apply(normalize_svtype)

    df_snv=df[df.SVTYPE_NORM=="SNV"].copy()
    df_sv=df[df.SVTYPE_NORM=="SV"].copy()

    print("SNV:",len(df_snv))
    print("SV :",len(df_sv))

    print("[6] spliting")

    # 1. Перемешиваем строки
    df_shuffled = df.sample(frac=1, random_state=42).reset_index(drop=True)

    # 2. Делим на train и temp (val+test), сохраняя пропорции классов
    train_df, temp = train_test_split(
        df_shuffled,
        test_size=0.7,            # 30% уйдёт на val+test
        stratify=df_shuffled["TRUE_VARIANT"],
        random_state=42
    )

    # 3. Делим temp на val и test, тоже со стратификацией
    val_df, test_df = train_test_split(
        temp,
        test_size=0.5,            # половина temp → test, половина → val
        stratify=temp["TRUE_VARIANT"],
        random_state=42
    )

    print(len(train_df), len(val_df), len(test_df))
    
    train_df[f"ML_PROB"] = 0.0
    test_df[f"ML_PROB"] = 0.0
    val_df[f"ML_PROB"] = 0.0

    features = [c for c in df.columns if c not in ["TRUE_VARIANT", "CHROM", "POS", "END", "REF", "ALT", "SVTYPE", "SVTYPE_NORM"]]

    print(f"\n=== Training model")

    X_train = train_df[features]
    y_train = train_df["TRUE_VARIANT"]
    X_val   = val_df[features]
    y_val   = val_df["TRUE_VARIANT"]
   
    # -----------------------------
    # 2. Обучение и валидация
    # -----------------------------
    results = []
    cat_cols = X_train.select_dtypes(include="object").columns.tolist() 

    param_grid = [
        dict(random_seed=254, l2_leaf_reg=5, depth=4, eval_metric="AUC", loss_function="Logloss"),
        dict(random_seed=487, l2_leaf_reg=0, depth=6, eval_metric="PRAUC", loss_function="Logloss"),
        dict(random_seed=789, l2_leaf_reg=2, depth=8, eval_metric="AUC", loss_function="Logloss", subsample=0.7, rsm=0.8, random_strength=1, bagging_temperature=1),
    ]

    for i, params in enumerate(param_grid):
        print(f"\n[{i+1}/{len(param_grid)}] Обучение модели с параметрами: {params}")

        model = CatBoostClassifier(
            iterations=100,
            early_stopping_rounds=20,
            learning_rate=0.03,
            verbose=100,
            posterior_sampling=True,
            **params
        )

        model.fit(X_train, y_train, eval_set=(X_val, y_val), cat_features=cat_cols, use_best_model=True)

        # вероятности на валидации
        val_prob = model.predict_proba(X_val)[:, 1]

        # подбор threshold под precision >= 0.9
        precision, recall, thresholds = precision_recall_curve(y_val, val_prob)

        target_precision = 0.99
        best_recall = 0
        best_threshold = 0

        for p, r, t in zip(precision, recall, np.append(thresholds, 1.0)):
            if p >= target_precision and r > best_recall:
                best_recall = r
                best_threshold = t

        print(f"target_precision = {target_precision}, best_recall = {best_recall}, best_threshold = {best_threshold}")
        fi = model.get_feature_importance(prettified=True)
        print(fi.head(8))

        results.append({
            "params": params,
            "model": model,
            "recall": best_recall,
            "threshold": best_threshold
        })
    
    save_hyperparam_results(results, f"{args.prefix}_hyperparams.tsv")

    best_result = max(results, key=lambda x: x["recall"])

    print("Лучшая модель:", best_result["params"])

    print(f"[7] prediction for model")

    test_prob_one = best_result["model"].predict_proba(test_df[features])[:, 1]

    test_pred_one = (test_prob_one >= best_result["threshold"]).astype(int)

    print("Отчёт для одной модели на test\n\n", classification_report(test_df["TRUE_VARIANT"], test_pred_one))

    df["ML_PROB"] = best_result["model"].predict_proba(df[features])[:, 1]
    df["ML_PRED"] = (df["ML_PROB"] >= best_result["threshold"]).astype(int)

    print("Финальный отчёт для одной модели на всём геноме\n\n", classification_report(df["TRUE_VARIANT"], df["ML_PRED"]))
   
    out_tsv = f"{args.prefix}.variants.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved {out_tsv}")

    # Save VCFs
    print(f"[8] VCF writing for ensemble")

    df_snv_sub = df[df.SVTYPE_NORM == "SNV"]

    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}.one_model.vcf", best_result["threshold"], "ML_PROB")

    best_result["model"].save_model(f"model_{args.prefix}.cbm")

    print("Model training and prediction complete.")

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

    y_pred_proba = best_result["model"].predict_proba(X_val)[:, 1]

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
    cm = confusion_matrix(df["TRUE_VARIANT"], df["ML_PRED"], labels=[0, 1])
    labels = ["Negative", "Positive"]

    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=labels, yticklabels=labels)
    plt.xlabel("Predicted label")
    plt.ylabel("True label")
    plt.title(f"Confusion Matrix (threshold = {best_result['threshold']})")
    plt.savefig(f"matrix_{args.prefix}.png", dpi=200)


if __name__=="__main__":

    p=argparse.ArgumentParser()

    p.add_argument("--mode", choices=["train", "predict"], default="predict")
    p.add_argument("--vcfs",nargs="+",required=True)
    p.add_argument("--truth_vcf")
    p.add_argument("--ml_model", help="Pre-trained CatBoost .cbm model for predict mode")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--merfin_pass_vcf")
    p.add_argument("--merfin_scores_bed")
    p.add_argument("--merqury")
    p.add_argument("--quast")
    p.add_argument("--repeat_gff")
    p.add_argument("--liftoff_gff")
    p.add_argument("--low_complex")
    p.add_argument("--flagger")
    p.add_argument("--out_vcf")
    p.add_argument("--out_table")
    p.add_argument("--exclude_telomers", help="BED file with telomere regions to exclude")
    p.add_argument("--prefix",required=True)


    args = p.parse_args()
    if args.mode == "predict" and not args.ml_model:
        p.error("--ml_model is required in predict mode and must point to a CatBoost .cbm file")
    main(args)
