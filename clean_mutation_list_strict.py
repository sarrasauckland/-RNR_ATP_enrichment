import pandas as pd
import re

# Load file
df = pd.read_csv("mutated_CDS_fixed.csv", header=None, names=["raw"])

# Regex pattern for valid gene names (E. coli style)
pattern = r"\b([a-z]{3,5}[A-Z]?[a-z]?\d*)\b"

# Extract possible matches
genes = []
for line in df["raw"].astype(str):
    found = re.findall(pattern, line)
    for g in found:
        # Filter out common false positives
        if g.lower() not in {
            "type","like","protein","domain","family","hypothetical",
            "putative","transporter","dependent","autotransporter",
            "prepilin","bifunctional","cds","system"
        } and len(g) <= 6:
            genes.append(g.lower())

# Deduplicate and save
unique = sorted(set(genes))
out = pd.DataFrame(unique, columns=["gene"])
out.to_csv("mutation_genes_clean_strict.csv", index=False)
print(f"✅ Cleaned {len(out)} genes saved to mutation_genes_clean_strict.csv")

