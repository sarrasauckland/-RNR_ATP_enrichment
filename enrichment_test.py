import pandas as pd
from scipy.stats import fisher_exact

# Define totals
total_genes = 4448
atp_genes = 296

# Load your SNP list
snps = pd.read_csv("mutated_CDS.csv")
snp_genes = snps["CDS"].str.extract(r"([A-Za-z0-9_]+)")[0].nunique()

# Load overlap file (SNPs found in ATP genes)
merged = pd.read_csv("snp_hits_in_atp_genes.csv")
overlap = merged["gene_name"].nunique()

# Build contingency table
# [[in_ATP_with_SNP, in_ATP_without_SNP],
#  [not_ATP_with_SNP, not_ATP_without_SNP]]
table = [
    [overlap, atp_genes - overlap],
    [snp_genes - overlap, total_genes - atp_genes - (snp_genes - overlap)]
]

oddsratio, p = fisher_exact(table, alternative="greater")

print("📊 Fisher’s Exact Test for ATP-gene enrichment")
print(f"Total genes: {total_genes}")
print(f"ATP genes: {atp_genes}")
print(f"Genes with SNPs: {snp_genes}")
print(f"ATP genes with SNPs: {overlap}")
print(f"\nContingency table: {table}")
print(f"Odds ratio: {oddsratio:.2f}")
print(f"P-value: {p:.3e}")

