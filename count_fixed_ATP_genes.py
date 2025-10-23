import pandas as pd

# --- Load your data ---
mutations = pd.read_csv("mutated_CDS_fixed.csv")  # has column 'CDS'
atp = pd.read_csv("atp_genes_with_names.csv")     # has columns 'b_number' and 'gene_name'

# --- Clean up gene name columns ---
mutations["CDS"] = mutations["CDS"].str.lower().str.strip()
atp["gene_name"] = atp["gene_name"].str.lower().str.strip()

# --- Find overlap between mutated genes and ATP-consuming genes ---
mutated_genes = set(mutations["CDS"])
atp_genes = set(atp["gene_name"])

overlap = mutated_genes.intersection(atp_genes)

# --- Output summary ---
print("🧬 FIXED SNP GENE-LEVEL COUNTS")
print("----------------------------------")
print(f"Total unique mutated genes: {len(mutated_genes)}")
print(f"ATP-consuming genes total: {len(atp_genes)}")
print(f"ATP-consuming genes mutated: {len(overlap)}")
print(f"Non-ATP genes mutated: {len(mutated_genes) - len(overlap)}")

# --- Optionally save overlapping genes to file ---
pd.DataFrame(sorted(overlap), columns=["gene_name"]).to_csv("mutated_ATP_fixed_overlap.csv", index=False)
print("\n✅ Saved overlapping ATP gene list to 'mutated_ATP_fixed_overlap.csv'")

