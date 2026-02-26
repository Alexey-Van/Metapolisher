def metrics(TP, FP, FN, TN):
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    accuracy = (TP + TN) / (TP + TN + FP + FN) if (TP + TN + FP + FN) > 0 else 0
    f1 = 2 * TP / (2 * TP + FP + FN) if (2 * TP + FP + FN) > 0 else 0
    fpr = FP / (FP + TN) if (FP + TN) > 0 else 0
    fnr = FN / (FN + TP) if (FN + TP) > 0 else 0

    return {
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "accuracy": accuracy,
        "f1": f1,
        "fpr": fpr,
        "fnr": fnr
    }

# пример
print(metrics(TP=231, FP=287, FN=83, TN=38792))
