import pandas as pd

# Load your files
snps = pd.read_csv("mutated_CDS.csv")            # ← your SNPs file
atp = pd.read_csv("atp_genes_with_names.csv")    # ← your ATP gene list

# Extract gene name from the CDS field
snps["gene_name"] = snps["CDS"].str.extract(r"([A-Za-z0-9_]+)")

# Merge on gene_name
merged = pd.merge(snps, atp, on="gene_name", how="inner")

print(f"✅ Found {len(merged)} SNPs in ATP-consuming genes.")
merged.to_csv("snp_hits_in_atp_genes.csv", index=False)
print("✅ Saved to snp_hits_in_atp_genes.csv")

