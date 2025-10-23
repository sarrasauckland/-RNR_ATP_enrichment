import pandas as pd
from scipy.stats import fisher_exact

# === CONSTANTS ===
TOTAL_GENES = 4448
ATP_TOTAL = 225

# === LOAD DATA ===
df = pd.read_csv("fixed_SNP_bylines_with_ATPflag.csv")
df.columns = ["line", "gene", "is_ATP_gene"]

# Clean formatting
df["line"] = df["line"].astype(str).str.strip().str.lower()
df["gene"] = df["gene"].astype(str).str.strip().str.lower()
df["is_ATP_gene"] = df["is_ATP_gene"].astype(str).str.lower().isin(["true", "t", "1"])

df = df.drop_duplicates(subset=["line", "gene"])

# === PER-LINE ENRICHMENT ===
results = []

for line, group in df.groupby("line"):
    total_mut = len(group)
    atp_mut = group["is_ATP_gene"].sum()

    # Contingency table
    a = atp_mut
    b = ATP_TOTAL - a
    c = total_mut - a
    d = (TOTAL_GENES - ATP_TOTAL) - c

    odds, p = fisher_exact([[a, b], [c, d]])

    results.append({
        "line": line,
        "total_mutations": total_mut,
        "atp_mutations": atp_mut,
        "odds_ratio": odds,
        "p_value": p
    })

# === OUTPUT ===
results_df = pd.DataFrame(results).sort_values("p_value", ascending=True)
results_df["significant"] = results_df["p_value"] < 0.05

results_df.to_csv("linewise_ATP_enrichment.csv", index=False)

print("📊 ATP Enrichment by Line")
print(results_df.round(4))
print("\n✅ Saved: linewise_ATP_enrichment.csv")

