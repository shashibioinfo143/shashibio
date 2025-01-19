import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the file
file_path = "C:/Users/shash/Downloads/PupilBioTest_PMP_revA.csv"
df = pd.read_csv(file_path)

# Display basic info and the first few rows to understand the structure

coverage_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
df['Total_CpG_Coverage'] = df[coverage_columns].sum(axis=1)

summary_stats = df.groupby("Tissue")["Total_CpG_Coverage"].agg(['median', 'mean', 'std'])
summary_stats["CV"] = summary_stats["std"] / summary_stats["mean"]

# Selecting only the median and CV for output
summary_stats = summary_stats[['median', 'CV']]
summary_stats

print(summary_stats)

# Set plot style
sns.set(style="whitegrid")

# Create a figure with subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Boxplot of CpG coverage per tissue
sns.boxplot(x="Tissue", y="Total_CpG_Coverage", data=df, ax=axes[0])
axes[0].set_title("CpG Coverage Distribution by Tissue")
axes[0].set_yscale("log")  # Log scale for better visibility

# Histogram of CpG coverage per tissue
sns.histplot(data=df, x="Total_CpG_Coverage", hue="Tissue", bins=50, log_scale=True, ax=axes[1], kde=True)
axes[1].set_title("Histogram of CpG Coverage by Tissue")

# Density plot of CpG coverage per tissue
sns.kdeplot(data=df, x="Total_CpG_Coverage", hue="Tissue", log_scale=True, ax=axes[2])
axes[2].set_title("Density Plot of CpG Coverage by Tissue")

# Adjust layout and display the plots
plt.tight_layout()
plt.show()
