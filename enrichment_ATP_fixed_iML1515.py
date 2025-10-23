from scipy.stats import fisher_exact

# --- Updated INPUT DATA ------------------------------------------------
# Model-based totals (iML1515 metabolic model)
total_genes = 1516          # total model genes
total_atp = 296             # ATP-consuming genes (unchanged)

# Fixed SNP summary from your previous analysis
mutated_genes = 860          # total genes with ≥1 fixed SNP
mutated_atp_genes = 76       # of those, how many are ATP-consuming
# -----------------------------------------------------------------------

# 1️⃣ Construct contingency table
# a = ATP genes mutated
# b = ATP genes NOT mutated
# c = non-ATP genes mutated
# d = non-ATP genes NOT mutated
a = mutated_atp_genes
b = total_atp - mutated_atp_genes
c = mutated_genes - mutated_atp_genes
d = total_genes - (a + b + c)

table = [[a, b],
         [c, d]]

# 2️⃣ Fisher’s exact test
odds_ratio, p_value = fisher_exact(table)

# 3️⃣ Output
print("GENE-LEVEL ENRICHMENT (ATP genes, fixed SNPs only, iML1515 model)")
print("==================================================================")
print(f"Total genes in model: {total_genes}")
print(f"ATP-consuming genes: {total_atp}")
print(f"Genes with ≥1 fixed SNP: {mutated_genes}")
print(f"ATP genes with fixed SNP: {mutated_atp_genes}")
print()
print(f"Contingency table = {table}")
print(f"Odds ratio = {odds_ratio:.2f}")
print(f"p-value = {p_value:.3e}")

