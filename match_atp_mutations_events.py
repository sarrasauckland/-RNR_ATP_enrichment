import pandas as pd

# Load event-level mutation list
mut = pd.read_csv("mutated_CDS.csv", names=["CDS"])
atp = pd.read_csv("atp_genes_with_names.csv", dtype=str)

# Clean and normalize
mut["gene_name"] = mut["CDS"].str.replace(" CDS", "", regex=False).str.strip()
atp["gene_name"] = atp["gene_name"].astype(str).str.strip()

# Event-level merge — keep duplicates (each line = one event)
merged = pd.merge(mut, atp, on="gene_name", how="inner")

# --- Event-level counts ---
total_events = len(mut)
atp_events = len(merged)
unique_genes = mut["gene_name"].nunique()
unique_atp_genes = atp["gene_name"].nunique()
unique_hit_genes = merged["gene_name"].nunique()

print("📊 Mutation Event Summary")
print(f"Total mutation events: {total_events}")
print(f"ATP-consuming mutation events: {atp_events}")
print(f"Unique mutated genes: {unique_genes}")
print(f"ATP-consuming genes total: {unique_atp_genes}")
print(f"ATP-consuming genes mutated: {unique_hit_genes}")

merged.to_csv("mutated_CDS_in_ATP_events.csv", index=False)
print("\n✅ Saved event-level matches to mutated_CDS_in_ATP_events.csv")

