import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(10, 6))

bins = np.arange(0, 301, 5)  # такие же бины как в статье (0–300, шаг 5)

plt.hist(v1, bins=bins, alpha=0.6, label="v1.0", color="#E69F00")
plt.hist(racon, bins=bins, alpha=0.6, label="Racon", color="#009E73")
plt.hist(racon_merfin, bins=bins, alpha=0.6, label="Racon + Merfin", color="#0072B2")

plt.yscale("log")  # логарифмическая ось Y

plt.xlabel("No. of corrections in 1-Mb windows", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.legend()

# Добавляем «Polishing hotspots»
plt.axvspan(50, 300, color="gray", alpha=0.15)
plt.text(160, 1.2, "Polishing hotspots", ha="center", fontsize=12)

plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Пример данных
df = pd.DataFrame({
    "CHROM": ["1","2","3","4","5","6","7","8","9","10",
              "11","12","13","14","15","16","17","18","19","20",
              "21","22","X","M"],
    "v09": [...],   # сюда твои данные
    "v10": [...],
})

# Устанавливаем порядок хромосом
df["CHROM"] = pd.Categorical(df["CHROM"],
    categories=[str(i) for i in range(1,23)] + ["X","M"],
    ordered=True
)
df = df.sort_values("CHROM")

x = np.arange(len(df))

plt.figure(figsize=(12,6))

plt.bar(x - 0.2, df["v09"], width=0.4, color="lightgray", label="v0.9")
plt.bar(x + 0.2, df["v10"], width=0.4, color="black", label="v1.0")

plt.xticks(x, df["CHROM"])
plt.xlabel("Chromosome")
plt.ylabel("Errors")
plt.legend()

plt.tight_layout()
plt.show()
