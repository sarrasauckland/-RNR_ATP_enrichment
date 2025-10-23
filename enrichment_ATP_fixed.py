import pandas as pd
from scipy.stats import fisher_exact

# --- Load datasets ---
mutations = pd.read_csv("mutated_CDS_fixed.csv")        # has 'CDS'
atp = pd.read_csv("atp_genes_with_names.csv")           # has 'gene_name'

# --- Clean and normalise names ---
def clean_name(x):
    if isinstance(x, str):
        x = x.lower().strip()
        x = x.replace(" cds", "")           # remove trailing CDS
        x = x.replace(" protein", "")       # remove generic labels
    return x

mutations["CDS"] = mutations["CDS"].apply(clean_name)
atp["gene_name"] = atp["gene_name"].apply(clean_name)

# --- Identify overlaps ---
mutated_genes = set(mutations["CDS"])
atp_genes = set(atp["gene_name"])
overlap = mutated_genes.intersection(atp_genes)

# --- Gene-level summary ---
total_genes = 4448   # total protein-coding genes in E. coli
total_atp = len(atp_genes)
mutated_gene_count = len(mutated_genes)
mutated_atp_count = len(overlap)

# --- Event-level summary ---
mutations["is_atp"] = mutations["CDS"].isin(atp_genes)
atp_event_count = mutations["is_atp"].sum()
total_events = len(mutations)
non_atp_events = total_events - atp_event_count

# --- Fisher’s exact tests ---
# Gene-level
a = mutated_atp_count
b = total_atp - mutated_atp_count
c = mutated_gene_count - mutated_atp_count
d = total_genes - (a + b + c)
odds_gene, p_gene = fisher_exact([[a, b], [c, d]])

# Event-level
non_atp_genes = total_genes - total_atp
odds_event, p_event = fisher_exact([[atp_event_count, non_atp_events],
                                   [total_atp, non_atp_genes]])

# --- Output summary ---
print("📊 ENRICHMENT SUMMARY (ATP-fixed SNPs)")
print("========================================")
print("EVENT-LEVEL DATA:")
print(f"  Total mutation events: {total_events}")
print(f"  ATP-consuming mutation events: {atp_event_count}")
print(f"  Non-ATP mutation events: {non_atp_events}")
print()
print("GENE-LEVEL DATA:")
print(f"  Total unique mutated genes: {mutated_gene_count}")
print(f"  ATP genes total: {total_atp}")
print(f"  ATP genes mutated: {mutated_atp_count}")
print(f"  Non-ATP genes mutated: {mutated_gene_count - mutated_atp_count}")
print()
print("📈 STATISTICAL RESULTS")
print(f"  Gene-level enrichment:  odds ratio = {odds_gene:.2f},  p = {p_gene:.3e}")
print(f"  Event-level enrichment: odds ratio = {odds_event:.2f}, p = {p_event:.3e}")

# --- Save overlap list ---
pd.DataFrame(sorted(overlap), columns=["gene_name"]).to_csv("mutated_ATP_fixed_overlap.csv", index=False)
print("\n✅ Saved overlapping ATP gene list to 'mutated_ATP_fixed_overlap.csv'")

