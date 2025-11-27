import pandas as pd
import matplotlib.pyplot as plt

#input the file

file_path = "C:/Users/shash/Downloads/PGS001298_hmPOS_GRCh38.txt"

# Loading and skip metadata lines starting with "#"
df = pd.read_csv(file_path, sep="\t", comment="#", dtype=str)

# Drop completely empty rows
df = df.dropna(how="all")

# Convert numeric columns
df["effect_weight"] = pd.to_numeric(df["effect_weight"], errors="coerce")
df["hm_chr"] = pd.to_numeric(df["hm_chr"], errors="coerce")
df["hm_pos"] = pd.to_numeric(df["hm_pos"], errors="coerce")

# -----------------------------
# Filter for chromosome 21 using hm_chr & hm_pos
# -----------------------------
df_chr21 = df[df["hm_chr"] == 21]

print("Number of variants on chromosome 21:", len(df_chr21))

# -----------------------------
# Plot histogram
# -----------------------------
plt.figure(figsize=(6,4))
plt.hist(df_chr21["effect_weight"].dropna(), bins=20)
plt.title("Distribution of effect_weight for Chromosome 21 (hm_chr/hm_pos)")
plt.xlabel("effect_weight")
plt.ylabel("Frequency")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()