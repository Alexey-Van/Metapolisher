from catboost import CatBoostClassifier
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

from catboost import CatBoostClassifier

help(CatBoostClassifier)


df = pd.read_csv("/mnt/data/vanichkinao/MetaPolisher/tmp/best_4__final_assembly.variants.tsv", sep="\t")
df = df.fillna(str(0.0))
features = [c for c in df.columns if c not in ["TRUE_VARIANT", "CHROM", "POS", "END", "REF", "ALT", "SVTYPE", "SVTYPE_NORM"]]

print(f"\n=== Training model")

X_train = df[features]
y_train = df["TRUE_VARIANT"]
X_val   = df[features]
y_val   = df["TRUE_VARIANT"]

model = CatBoostClassifier(
            iterations=100,
            learning_rate=0.03,
            early_stopping_rounds=20,
            # subsample=0.7,
            # rsm=0.8,
            # random_strength=2,
            # bagging_temperature=2,
            verbose=100,
            posterior_sampling=True,
            eval_metric="AUC", 
            loss_function="Logloss"
        )

model.fit(X_train, y_train, eval_set=(X_val, y_val), use_best_model=True)