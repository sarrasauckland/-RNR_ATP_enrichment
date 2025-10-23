import pandas as pd
from scipy.stats import fisher_exact

TOTAL_GENES = 4448
ATP_TOTAL = 225
FILE = "fixed_SNP_bylines_with_ATPnew.csv"

df = pd.read_csv(FILE, header=None, names=["line", "gene", "is_ATP_gene"])
df["line"] = df["line"].astype(str).str.strip().str.lower().str.replace(" ", "")
df["gene"] = df["gene"].astype(str).str.strip().str.lower()
df["is_ATP_gene"] = df["is_ATP_gene"].astype(str).str.lower().isin(["true","t","1"])

print("Detected lines:", sorted(df["line"].unique()), "\n")

for line, group in df.groupby("line", sort=True):
    total_events = len(group)
    atp_events = int(group["is_ATP_gene"].sum())
    non_atp_events = total_events - atp_events

    unique_genes = group.drop_duplicates(subset=["gene"])
    a = int(unique_genes["is_ATP_gene"].sum())
    c = len(unique_genes) - a
    b = ATP_TOTAL - a
    d = (TOTAL_GENES - ATP_TOTAL) - c

    odds_gene, p_gene = fisher_exact([[a,b],[c,d]])
    odds_event, p_event = fisher_exact([[atp_events, non_atp_events],
                                        [ATP_TOTAL, TOTAL_GENES - ATP_TOTAL]])

    print(f"🧩 {line.upper()}")
    print("========================================")
    print("EVENT-LEVEL DATA:")
    print(f"  Total mutation events: {total_events}")
    print(f"  ATP-consuming mutation events: {atp_events}")
    print(f"  Non-ATP mutation events: {non_atp_events}\n")

    print("GENE-LEVEL DATA:")
    print(f"  Total unique mutated genes: {len(unique_genes)}")
    print(f"  ATP genes total: {ATP_TOTAL}")
    print(f"  ATP genes mutated: {a}")
    print(f"  Non-ATP genes mutated: {c}\n")

    print("📈 STATISTICAL RESULTS")
    print(f"  Gene-level enrichment:  odds ratio = {odds_gene:.2f},  p = {p_gene:.3e}")
    print(f"  Event-level enrichment: odds ratio = {odds_event:.2f},  p = {p_event:.3e}")
    print("========================================\n")

total_all = len(df)
atp_all = int(df["is_ATP_gene"].sum())
non_atp_all = total_all - atp_all
unique_all = df.drop_duplicates(subset=["gene"])
a = int(unique_all["is_ATP_gene"].sum())
c = len(unique_all) - a
b = ATP_TOTAL - a
d = (TOTAL_GENES - ATP_TOTAL) - c
odds_gene_all, p_gene_all = fisher_exact([[a,b],[c,d]])
odds_event_all, p_event_all = fisher_exact([[atp_all, non_atp_all],
                                            [ATP_TOTAL, TOTAL_GENES - ATP_TOTAL]])

print("🌍 ALL LINES COMBINED")
print("========================================")
print("EVENT-LEVEL DATA:")
print(f"  Total mutation events: {total_all}")
print(f"  ATP-consuming mutation events: {atp_all}")
print(f"  Non-ATP mutation events: {non_atp_all}\n")

print("GENE-LEVEL DATA:")
print(f"  Total unique mutated genes: {len(unique_all)}")
print(f"  ATP genes total: {ATP_TOTAL}")
print(f"  ATP genes mutated: {a}")
print(f"  Non-ATP genes mutated: {c}\n")

print("📈 STATISTICAL RESULTS")
print(f"  Gene-level enrichment:  odds ratio = {odds_gene_all:.2f},  p = {p_gene_all:.3e}")
print(f"  Event-level enrichment: odds ratio = {odds_event_all:.2f},  p = {p_event_all:.3e}")
print("========================================")

