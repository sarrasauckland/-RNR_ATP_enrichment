from scipy.stats import fisher_exact
import pandas as pd

# === Load and clean ===
df = pd.read_csv("fixed_SNP_bylines_with_ATPnew.csv", header=None, skiprows=1,
                 names=["line", "gene", "is_ATP_gene"])
df["line"] = (
    df["line"].astype(str)
    .str.strip()
    .str.lower()
    .replace(r"\s+", " ", regex=True)
    .replace("line ", "line_", regex=True)
)
df["gene"] = df["gene"].astype(str).str.strip().str.lower()
df["is_ATP_gene"] = df["is_ATP_gene"].astype(str).str.lower().isin(["true", "t", "1"])

ATP_TOTAL = 225
GENOME_TOTAL = 4448

results = []

# === Loop over lines ===
for line, subset in df.groupby("line"):
    total_mut = len(subset)
    atp_mut = subset["is_ATP_gene"].sum()
    non_atp_mut = total_mut - atp_mut

    a = atp_mut
    b = ATP_TOTAL - atp_mut
    c = non_atp_mut
    d = (GENOME_TOTAL - ATP_TOTAL) - non_atp_mut

    odds, p = fisher_exact([[a, b], [c, d]])

    results.append({
        "line": line,
        "total_mutations": total_mut,
        "atp_mutations": atp_mut,
        "odds_ratio": odds,
        "p_value": p
    })

# === Output ===
results_df = pd.DataFrame(results).sort_values("p_value", ascending=True)
print("\n📊 ATP ENRICHMENT BY LINE")
print(results_df.round(4))
results_df.to_csv("linewise_ATP_enrichment_fixed.csv", index=False)
print("\n✅ Saved: linewise_ATP_enrichment_fixed.csv")

