import pandas as pd
import re

# Load your mutation list
df = pd.read_csv("mutated_CDS_fixed.csv", header=None, names=["raw"])

# Extract likely gene names (matches typical E. coli gene formats: e.g. cyaA, yfgL, iscS)
df["gene"] = df["raw"].str.extract(r"\b([a-z]{2,5}[A-Z]?[a-z]*\d*)\b")

# Drop blanks and clean up
df = df.dropna(subset=["gene"])
df["gene"] = df["gene"].str.lower().str.strip()
df_clean = df[["gene"]].drop_duplicates()

# Save cleaned list
df_clean.to_csv("mutation_genes_clean.csv", index=False)
print(f"✅ Cleaned {len(df_clean)} gene entries saved to mutation_genes_clean.csv")

