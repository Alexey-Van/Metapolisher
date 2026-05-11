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
            if abs(pos[i] - t[j]) <= 1:
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

    exclude_trees = None
    if args.exclude_telomers:
        exclude_trees = load_telomers(args.exclude_telomers)


    df=outer_merge_vcfs(args.vcfs, exclude_trees)

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
        test_size=0.3,            # 30% уйдёт на val+test
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

    # chroms_present = [c for c in df["CHROM"].unique()]

    # train_chroms = chroms_present[:10]      # chr1–chr16
    # val_chroms   = chroms_present[10:17]    # chr17–chr19
    # test_chroms  = chroms_present[17:]      # chr20–chr22, X, Y, M

    # train_df = df[df["CHROM"].isin(train_chroms)]
    # val_df   = df[df["CHROM"].isin(val_chroms)]
    # test_df  = df[df["CHROM"].isin(test_chroms)]

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
    models = []

    cat_cols = X_train.select_dtypes(include="object").columns.tolist() 

    param_grid = [
        dict(random_seed=11),
        dict(random_seed=22),
        dict(random_seed=33),
        dict(random_seed=44),
        dict(random_seed=55),
        dict(random_seed=23),
        dict(random_seed=42),
        dict(random_seed=43),
        dict(random_seed=58),
        dict(random_seed=102),
    ]

    for i, params in enumerate(param_grid):
        print(f"\n[{i+1}/{len(param_grid)}] Обучение модели с параметрами: {params}")

        model = CatBoostClassifier(
            iterations=300,
            learning_rate=0.03,
            loss_function="Logloss",
            l2_leaf_reg=10,
            auto_class_weights="SqrtBalanced",
            depth=6,
            subsample=0.7,
            rsm=0.8,
            random_strength=2,
            bagging_temperature=2,
            eval_metric="PRAUC",
            verbose=100,
            posterior_sampling=True,
            **params
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

        print(f"target_precision = {target_precision}, best_recall = {best_recall}, best_threshold = {best_threshold}")

        results.append({
            "params": params,
            "model": model,
            "recall": best_recall,
            "threshold": best_threshold
        })
    
    save_hyperparam_results(results, f"{args.prefix}_hyperparams.tsv")

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

    # извлекаем модели и их recall
    models = [x["model"] for x in top_models]
    recalls = np.array([x["recall"] for x in top_models])

    # нормируем веса
    weights = recalls / recalls.sum()

    print("Weights:", weights)


    test_prob_one = top_models[1]["model"].predict_proba(test_df[features])[:, 1]

    test_pred_one = (test_prob_one >= top_models[1]["threshold"]).astype(int)

    print("Отчёт для одной модели\n\n", classification_report(test_df["TRUE_VARIANT"], test_pred_one))

    # вероятности всех моделей на test
    test_probs_list = [
        m.predict_proba(test_df[features])[:, 1]
        for m in models
    ]

    # взвешенное усреднение
    test_prob = np.average(test_probs_list, axis=0, weights=weights)

    # бинарные предсказания
    test_pred = (test_prob >= best_threshold).astype(int)

    print("Отчёт для ансамбля\n\n", classification_report(test_df["TRUE_VARIANT"], test_pred))

    # вероятности ансамбля по всему геному
    all_probs_list = [
        m["model"].predict_proba(df[features])[:, 1]
        for m in top_models
    ]

    probs = np.vstack(all_probs_list)          # (n_models, n_variants)
    w = np.array(weights, dtype=float)        # (n_models,)

    df["ML_PROB"] = np.average(probs, axis=0, weights=w)
    df["ML_PRED"] = (df["ML_PROB"] >= best_threshold).astype(int)


    print("Финальный отчёт для финального ансамбля\n\n", classification_report(df["TRUE_VARIANT"], df["ML_PRED"]))

    df.loc[df.index, f"ML_PROB_one"] = predict_block(df, top_models[1]["model"])
    df["ML_PRED_one"] = (df["ML_PROB_one"] >= top_models[1]["best_threshold"]).astype(int)

    print("Финальный отчёт для финального ансамбля\n\n", classification_report(df["TRUE_VARIANT"], df["ML_PRED_one"]))


   
    out_tsv = f"{args.prefix}.variants.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved {out_tsv}")

    # Save VCFs
    print(f"[8] VCF writing for ensemble")

    df_snv_sub = df[df.SVTYPE_NORM == "SNV"]

    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}.snv.vcf", best_threshold, ML_PROB)

    write_vcf(df_snv_sub, args.vcfs, f"{args.prefix}.snv.vcf", top_models[1]["best_threshold"], ML_PROB_one)

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