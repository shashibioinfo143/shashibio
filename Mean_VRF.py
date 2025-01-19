import pandas as pd

# Load the file
file_path = "C:/Users/shash/Downloads/PupilBioTest_PMP_revA.csv"
df = pd.read_csv(file_path)

coverage_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']

# Calculate total CpG coverage per row
df["Total_CpG_Coverage"] = df[coverage_columns].sum(axis=1)

# Compute VRF for each PMP (PMP read count / total coverage)
for col in coverage_columns:
    df[f"VRF_{col}"] = df[col] / df["Total_CpG_Coverage"]

# Reshape data to long format for grouping by PMP and Tissue
vrf_cols = [f"VRF_{col}" for col in coverage_columns]
df_vrf = df.melt(id_vars=["CpG_Coordinates", "Tissue"], value_vars=vrf_cols,
                 var_name="PMP", value_name="VRF")

# Clean PMP names to remove "VRF_" prefix
#df_vrf["PMP"] = df_vrf["PMP"].str.replace("VRF_", "")

# Group by PMP and Tissue to compute mean VRF
mean_vrf = df_vrf.groupby(["PMP", "Tissue"])["VRF"].mean().unstack()

mean_vrf
print(mean_vrf)
