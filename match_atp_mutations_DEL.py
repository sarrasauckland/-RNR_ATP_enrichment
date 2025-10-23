import pandas as pd

# --- Load files ---
mut = pd.read_csv("mutated_CDS_DEL.csv", names=["CDS"])
atp = pd.read_csv("atp_genes_with_names.csv", dtype=str)

# --- Clean names ---
mut["gene_name"] = mut["CDS"].str.replace(" CDS", "", regex=False).str.strip()
atp["gene_name"] = atp["gene_name"].astype(str).str.strip()

# --- Match events (keep duplicates = count every mutation) ---
merged = pd.merge(mut, atp, on="gene_name", how="inner")

# --- Count stats ---
total_events = len(mut)
atp_events = len(merged)
unique_mut_genes = mut["gene_name"].nunique()
unique_atp_genes = atp["gene_name"].nunique()
unique_hit_genes = merged["gene_name"].nunique()

# --- Output ---
print("📊 Mutation Event Summary (DEL version)")
print(f"Total mutation events: {total_events}")
print(f"ATP-consuming mutation events: {atp_events}")
print(f"Unique mutated genes: {unique_mut_genes}")
print(f"ATP-consuming genes total: {unique_atp_genes}")
print(f"ATP-consuming genes mutated: {unique_hit_genes}")

merged.to_csv("mutated_CDS_DEL_in_ATP.csv", index=False)
print("\n✅ Saved event-level matches to mutated_CDS_DEL_in_ATP.csv")

