# Dependencies
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the actual file
data = pd.read_csv("simrel_scores_output.tsv", sep="\t")

# Debug: Print columns to ensure they are read correctly
print("DataFrame columns:", data.columns)

# Strip leading/trailing whitespace from column names (just in case)
data.columns = data.columns.str.strip()

# Replace 'NA' with NaN and convert score columns to numeric
for col in ["MF_Score", "CC_Score", "BP_Score", "Total_Score"]:
    if col in data.columns:  # Ensure the column exists
        data[col] = pd.to_numeric(data[col].replace("NA", np.nan), errors='coerce')
    else:
        print(f"Warning: Column {col} is missing from the data!")

# Prepare data for boxplot
boxplot_data = [data[col].dropna() for col in ["MF_Score", "CC_Score", "BP_Score", "Total_Score"] if col in data.columns]

# Generate the boxplot
plt.figure(figsize=(10, 6))
plt.boxplot(
    boxplot_data,
    tick_labels=["MF_Score", "CC_Score", "BP_Score", "Total_Score"],  # Updated to tick_labels
    showmeans=True
)
plt.title("Boxplot of Scores")
plt.ylabel("Score")
plt.xlabel("Score Type")
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Save the plot as an image
plt.savefig("boxplot_scores.png", dpi=300)
print("Boxplot saved as 'boxplot_scores.png'.")
