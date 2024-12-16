import pandas as pd

# Read the  results file
input_file = "results.tsv"
output_file = "go_terms_DeepFRI.tsv"

# Load the file into a DataFrame
df = pd.read_csv(input_file, sep="\t", header=0, names=["Protein", "GO_term", "Score", "Annotation", "Neural_net", "DeepFRI_mode", "DB_hit", "DB_name", "Identity", "Coverage"])

# Group by Protein and aggregate GO terms
grouped = df.groupby("Protein")["GO_term"].apply(lambda terms: ",".join(terms.unique())).reset_index()

# Rename columns
grouped.columns = ["Protein", "GO_terms"]

# Save the output file
grouped.to_csv(output_file, sep="\t", index=False)

print(f"File saved as {output_file}")
