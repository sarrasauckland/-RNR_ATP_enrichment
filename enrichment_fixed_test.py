from scipy.stats import fisher_exact

# --- Observed counts (FIXED SNPs version) ---
total_genes = 4448          # all genes in REL606
total_atp = 296             # ATP-consuming genes
mutated_genes = 1025        # unique mutated genes (unchanged)
mutated_atp_genes = 78      # ATP-consuming genes mutated (unchanged unless new)

# --- Event-level counts ---
total_events = 1265          # total fixed mutation events
atp_events = 107             # change this if your new count differs

# -------------------------------
# GENE-LEVEL TEST
# -------------------------------
a = mutated_atp_genes                  # ATP genes mutated
b = total_atp - mutated_atp_genes      # ATP genes not mutated
c = mutated_genes - mutated_atp_genes  # Non-ATP genes mutated
d = total_genes - (a + b + c)
odds_gene, p_gene = fisher_exact([[a, b], [c, d]])

# -------------------------------
# EVENT-LEVEL TEST
# -------------------------------
mut_non_atp_events = total_events - atp_events
non_atp_genes = total_genes - total_atp
odds_event, p_event = fisher_exact([[atp_events, mut_non_atp_events],
                                    [total_atp, non_atp_genes]])

# -------------------------------
# OUTPUT
# -------------------------------
print("📊 Enrichment Summary (FIXED SNP version)")
print("----------------------------------------")
print("EVENT-LEVEL DATA:")
print(f"  Total mutation events: {total_events}")
print(f"  ATP-consuming mutation events: {atp_events}")
print(f"  Non-ATP mutation events: {mut_non_atp_events}")
print()
print("GENE-LEVEL DATA:")
print(f"  Total unique mutated genes: {mutated_genes}")
print(f"  ATP genes total: {total_atp}")
print(f"  ATP genes mutated: {mutated_atp_genes}")
print(f"  Non-ATP genes mutated: {mutated_genes - mutated_atp_genes}")
print()
print("📈 STATISTICAL RESULTS")
print(f"  Gene-level enrichment:  odds ratio = {odds_gene:.2f},  p = {p_gene:.3e}")
print(f"  Event-level enrichment: odds ratio = {odds_event:.2f}, p = {p_event:.3e}")

