import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------
input_file = "C:/Users/shash/Downloads/PGS001298_hmPOS_GRCh38.txt"
PLOT_DIR = "./plots/"
# ---------------------------

# Create directories if needed
import os
os.makedirs(PLOT_DIR, exist_ok=True)

# ---------------------------
# LOAD & BASIC CLEANING
# ---------------------------
df = pd.read_csv(input_file, sep="\t", comment="#", dtype=str)
df = df.dropna(how="all")

# Fill missing rsID using hm_rsID
if "hm_rsID" in df.columns:
    df["rsID"] = df["rsID"].fillna(df["hm_rsID"])

# Convert numeric columns
num_cols = ["chr_position", "effect_weight", "hm_pos"]
for col in num_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

print("\n===== BASIC INFO =====")
print(df.info())
print("\n===== HEAD =====")
print(df.head())

# ---------------------------
# SUMMARY STATISTICS
# ---------------------------
print("\n===== SUMMARY STATISTICS =====")
print(df.describe(include="all"))

# ---------------------------
# MISSINGNESS ANALYSIS
# ---------------------------
missing = df.isnull().mean().sort_values(ascending=False)
print("\n===== MISSINGNESS =====")
print(missing)

plt.figure(figsize=(8,5))
missing.plot(kind="bar")
plt.title("Missingness per Column")
plt.ylabel("Fraction Missing")
plt.tight_layout()
plt.savefig(PLOT_DIR + "missingness.png")
plt.close()

# ---------------------------
# EFFECT WEIGHT DISTRIBUTION
# ---------------------------
plt.figure(figsize=(8,5))
sns.histplot(df["effect_weight"].dropna(), bins=50, kde=True)
plt.title("Distribution of Effect Weights")
plt.xlabel("Effect Weight")
plt.tight_layout()
plt.savefig(PLOT_DIR + "effect_weight_distribution.png")
plt.close()

# ---------------------------
# CHROMOSOME-SPECIFIC ANALYSIS
# ---------------------------
if "chr_name" in df.columns:
    chr_counts = df["chr_name"].value_counts().sort_index()
    print("\n===== VARIANTS PER CHROMOSOME =====")
    print(chr_counts)

    plt.figure(figsize=(8,5))
    chr_counts.plot(kind="bar")
    plt.title("Number of Variants per Chromosome")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(PLOT_DIR + "variants_per_chromosome.png")
    plt.close()

# ---------------------------
# EFFECT ALLELE SUMMARY
# ---------------------------
if "effect_allele" in df.columns:
    effect_allele_counts = df["effect_allele"].value_counts()
    print("\n===== EFFECT ALLELE COUNTS =====")
    print(effect_allele_counts)

    plt.figure(figsize=(8,5))
    effect_allele_counts.plot(kind="bar")
    plt.title("Effect Allele Frequency")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(PLOT_DIR + "effect_allele_frequency.png")
    plt.close()

# ---------------------------
# OUTLIER ANALYSIS
# ---------------------------
q1 = df["effect_weight"].quantile(0.25)
q3 = df["effect_weight"].quantile(0.75)
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr

outliers = df[(df["effect_weight"] < lower) | (df["effect_weight"] > upper)]

print("\n===== OUTLIER SUMMARY =====")
print("IQR Lower Bound:", lower)
print("IQR Upper Bound:", upper)
print("Number of Outliers:", len(outliers))

# ---------------------------
# CORRELATION ANALYSIS
# ---------------------------
num_df = df[num_cols].dropna()

if len(num_df) > 0:
    plt.figure(figsize=(6,4))
    sns.heatmap(num_df.corr(), annot=True, cmap="coolwarm")
    plt.title("Correlation Matrix (Numeric Variables)")
    plt.tight_layout()
    plt.savefig(PLOT_DIR + "correlation_matrix.png")
    plt.close()

print("\nEDA Completed. Plots saved in:", PLOT_DIR)
