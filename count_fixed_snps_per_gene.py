import pandas as pd

# --- Load your fixed SNP file ---
df = pd.read_csv("mutated_CDS_fixed.csv")

# --- Clean up gene names ---
df["gene"] = df["CDS"].str.replace(" CDS", "", regex=False).str.strip().str.lower()

# --- Count SNPs per gene ---
snp_counts = df["gene"].value_counts().reset_index()
snp_counts.columns = ["gene", "snp_count"]

# --- Save the counts to a new CSV ---
snp_counts.to_csv("fixed_snp_counts.csv", index=False)

# --- Print summary ---
print(f"✅ Processed {len(df)} SNP entries across {len(snp_counts)} unique genes.")
print("Top 10 genes with most fixed SNPs:")
print(snp_counts.head(10))

