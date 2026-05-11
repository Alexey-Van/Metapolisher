# извлекаем модели и их recall
models = [x["model"] for x in top_models]
recalls = np.array([x["recall"] for x in top_models])

# нормируем веса
weights = recalls / recalls.sum()

print("Weights:", weights)

# вероятности всех моделей на test
test_probs_list = [
    m.predict_proba(test_df[features])[:, 1]
    for m in models
]

# взвешенное усреднение
test_prob = np.average(test_probs_list, axis=0, weights=weights)

# бинарные предсказания
test_pred = (test_prob >= best_threshold).astype(int)

print(classification_report(test_df["TRUE_VARIANT"], test_pred))
