import pandas as pd

# Load mutation events and ATP gene reference
mut = pd.read_csv("mutated_CDS.csv", names=["CDS"])
atp = pd.read_csv("atp_genes_with_names.csv", dtype=str)

# Clean text
mut["gene_name"] = mut["CDS"].str.replace(" CDS", "", regex=False).str.strip()
atp["gene_name"] = atp["gene_name"].astype(str).str.strip()

# Find matches (each mutated gene that’s ATP-consuming)
merged = pd.merge(mut, atp, on="gene_name", how="inner")

print(f"✅ Found {len(merged)} mutation events in ATP-consuming genes.")
merged.to_csv("mutated_CDS_in_ATP.csv", index=False)
print("✅ Saved to mutated_CDS_in_ATP.csv")

# Summary for enrichment
unique_mut_genes = mut['gene_name'].nunique()
unique_atp_genes = atp['gene_name'].nunique()
unique_hits = merged['gene_name'].nunique()

print(f"\nUnique mutated genes: {unique_mut_genes}")
print(f"ATP-consuming genes total: {unique_atp_genes}")
print(f"ATP-consuming genes mutated: {unique_hits}")

