import matplotlib.pyplot as plt
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

def main(args):

    df = pd.read_csv(args.tsv, sep="\t")

    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)


    plt.hist(df.loc[df.TRUE_VARIANT == 0, "ML_PROB"], bins=50, alpha=0.5, label="0")
    plt.hist(df.loc[df.TRUE_VARIANT == 1, "ML_PROB"], bins=50, alpha=0.5, label="1")
    plt.legend()
    plt.yscale("log")

    plt.savefig("hist.png", dpi=200)

    from sklearn.metrics import precision_recall_curve

    y_true = df.TRUE_VARIANT.values
    y_score = df.ML_PROB.values

    prec, rec, thr = precision_recall_curve(y_true, y_score)

    plt.plot(rec, prec)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.savefig("pr_curve.png", dpi=200)

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True, help="Input TSV from main pipeline")
    p.add_argument("--threshold",type=float,default=0.9)
    main(p.parse_args())