import pandas as pd

# ---------------------------
input_file = "C:/Users/shash/Downloads/PGS001298_hmPOS_GRCh38.txt"
output_file = "PGS001298_hmPOS_GRCh38_cleaned.csv"
# ----------------------------

def preprocess_pgs_file(input_file, output_file=None):
    # Load file, skip metadata lines beginning with "#"
    df = pd.read_csv(input_file, sep="\t", comment="#", dtype=str)

    # 1. Drop fully empty rows
    df = df.dropna(how="all")

    # 2. Fill missing rsID values with hm_rsID
    if "rsID" in df.columns and "hm_rsID" in df.columns:
        df["rsID"] = df["rsID"].fillna(df["hm_rsID"])

    # 3. Convert numeric columns if present
    numeric_cols = ["chr_position", "effect_weight", "hm_pos"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # 4. Remove rows missing essential identifiers
    df = df[
        df["rsID"].notna() |
        (df["chr_name"].notna() & df["chr_position"].notna())
    ]

    # 5. Reset index
    df = df.reset_index(drop=True)

    # 6. Save file
    if output_file:
        df.to_csv(output_file, index=False)

    return df


# ---------- Run preprocessing ----------
clean_df = preprocess_pgs_file(input_file, output_file)

print("Preprocessing complete!")
print("Cleaned shape:", clean_df.shape)
print(clean_df.head())
