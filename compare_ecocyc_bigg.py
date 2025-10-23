import pandas as pd

# Load both files
ecocyc = pd.read_csv("gene_counts_atp.tsv", sep="\t")
bigg = pd.read_csv("atp_genes_with_pathways.csv", on_bad_lines="skip")

# Normalize gene names
ecocyc["gene"] = ecocyc["gene"].str.lower().str.strip()
bigg["gene_name"] = bigg["gene_name"].str.lower().str.strip()

# Unique sets
ecocyc_genes = set(ecocyc["gene"].dropna())
bigg_genes = set(bigg["gene_name"].dropna())

# Compare
shared = ecocyc_genes & bigg_genes
only_ecocyc = ecocyc_genes - bigg_genes
only_bigg = bigg_genes - ecocyc_genes

# Summary
print("EcoCyc total genes:", len(ecocyc_genes))
print("BiGG total genes:", len(bigg_genes))
print("Shared genes:", len(shared))
print("Unique to EcoCyc:", len(only_ecocyc))
print("Unique to BiGG:", len(only_bigg))

# Save lists
pd.Series(sorted(only_ecocyc)).to_csv("unique_to_ecocyc.csv", index=False)
pd.Series(sorted(only_bigg)).to_csv("unique_to_bigg.csv", index=False)
pd.Series(sorted(shared)).to_csv("shared_genes.csv", index=False)

