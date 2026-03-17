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


# ------------------------------------------------------------
# utilities
# ------------------------------------------------------------

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
            "SVTYPE": info.get("SVTYPE","SNV"),
            "SVLEN": normalize_info_value("SVLEN", info.get("SVLEN")),
            f"{prefix}_present":1
        }

        for k,val in info.items():
            row[f"{prefix}_INFO_{k}"] = normalize_info_value(k,val)

        row.update(parse_format_fields(v,prefix))

        rows.append(row)

    return pd.DataFrame(rows)


def outer_merge_vcfs(paths):

    dfs = [load_vcf_as_df(p) for p in paths]

    df = dfs[0]

    for d in dfs[1:]:

        df = df.merge(
            d,
            on=["CHROM","POS","END","REF","ALT","SVTYPE","SVLEN"],
            how="outer"
        )

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

    truth={}

    for v in VCF(vcf):

        truth.setdefault(v.CHROM,[]).append(v.POS)

    for k in truth:
        truth[k]=np.array(truth[k])

    return truth


def build_truth_labels(df,truth):

    labels=np.zeros(len(df))

    for chrom,idx in df.groupby("CHROM").groups.items():

        if chrom not in truth:
            continue

        pos=df.loc[idx,"POS"].values
        t=truth[chrom]

        diff=np.abs(pos[:,None]-t)

        labels[idx]=(diff<=1).any(axis=1)

    return labels.astype(int)


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

    cat_cols=X.select_dtypes(include="object").columns.tolist()

    X_tr,X_te,y_tr,y_te=train_test_split(
        X,y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    model=CatBoostClassifier(
        iterations=300,
        depth=12,
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

    return model


def predict_block(df,model):

    drop_cols={
        "TRUE_VARIANT","CHROM","POS","END",
        "REF","ALT","SVTYPE","SVTYPE_NORM"
    }
    
    fi = model.get_feature_importance(prettified=True)
    print(fi.head(10))

    X=df.drop(columns=[c for c in drop_cols if c in df.columns])

    return model.predict_proba(X)[:,1]


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


# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main(args):

    print("[1] loading VCFs")

    df=outer_merge_vcfs(args.vcfs)

    # print("[2] merfin")

    # df=annotate_merfin_support(df,args.merfin_pass_vcf)

    # if args.merfin_scores_bed:
    #     df=annotate_merfin_scores(df,args.merfin_scores_bed)

    print("[3] annotations")

    repeat=load_gff3(args.repeat_gff,"Target")
    liftoff=load_gff3(args.liftoff_gff)
    low_complex=load_bed(args.low_complex)
    flagger=load_bed(args.flagger,3)

    annotate_interval(df,repeat,"repeat")
    annotate_interval(df,liftoff,"liftoff")
    annotate_interval(df,low_complex,"low_complex")
    annotate_interval(df,flagger,"flagger")

    print("[4] truth")

    truth=load_truth(args.truth_vcf)

    df["TRUE_VARIANT"]=build_truth_labels(df,truth)

    print("[5] SVTYPE")

    df["SVTYPE_NORM"]=df["SVTYPE"].apply(normalize_svtype)

    df_snv=df[df.SVTYPE_NORM=="SNV"].copy()
    df_sv=df[df.SVTYPE_NORM=="SV"].copy()

    print("SNV:",len(df_snv))
    print("SV :",len(df_sv))

    print("[6] training")

    model_snv=train_catboost(df_snv)
    model_sv=train_catboost(df_sv)

    print("[7] prediction")

    df["ML_PROB"] = 0.0

    df.loc[df_snv.index, "ML_PROB"] = predict_block(df_snv, model_snv)
    df.loc[df_sv.index,  "ML_PROB"] = predict_block(df_sv,  model_sv)

    df.to_csv(args.out_table, sep="\t", index=False)

    print("[8] VCF")

    df_snv = df[df.SVTYPE_NORM == "SNV"]
    df_sv  = df[df.SVTYPE_NORM == "SV"]

    write_vcf(df_snv, args.vcfs, "snv.filtered.vcf", args.threshold)
    write_vcf(df_sv,  args.vcfs, "sv.filtered.vcf",  args.threshold)
    write_vcf(df,     args.vcfs, args.out_vcf,       args.threshold)


if __name__=="__main__":

    p=argparse.ArgumentParser()

    p.add_argument("--vcfs",nargs="+",required=True)
    p.add_argument("--truth_vcf",required=True)
    p.add_argument("--merfin_pass_vcf")
    p.add_argument("--merfin_scores_bed")
    p.add_argument("--repeat_gff",required=True)
    p.add_argument("--liftoff_gff",required=True)
    p.add_argument("--low_complex",required=True)
    p.add_argument("--flagger",required=True)
    p.add_argument("--out_vcf",required=True)
    p.add_argument("--out_table",default="variant_features.tsv")
    p.add_argument("--threshold",type=float,default=0.9)

    main(p.parse_args())