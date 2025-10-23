import pandas as pd, re

# Input files
MAP = "reaction_gene_map_atp.tsv"
GENES = "genes_to_lookup.txt"

# Load the reaction map
df = pd.read_csv(MAP, sep="\t", dtype=str).fillna("")
df.columns = [c.strip().lower() for c in df.columns]

# Keep reactions where ATP is on the left (ATP consumed)
def left_has_atp(s: str) -> bool:
    parts = re.split(r'\s*(?:<[-−]*>|[-−]*>|→|⇌|<->|->)\s*', str(s), maxsplit=1)
    left = parts[0].upper() if parts else ""
    return bool(re.search(r'\bATP\b', left))

df = df[df["reaction"].apply(left_has_atp)].copy()

# Load gene list
with open(GENES) as f:
    targets = {ln.strip().lower() for ln in f if ln.strip()}

df["gene"] = df["gene"].astype(str).str.strip()
df = df[df["gene"].str.lower().isin(targets)].copy()

# Helper functions
def uniq(vals): 
    return sorted(v for v in set(vals) if v)

# Aggregate reactions per gene
out = (df.groupby("gene", as_index=False)
          .agg(
              n_reactions=("reaction", lambda s: len(set(s))),
              enzymes=("enzyme", lambda s: "; ".join(uniq(s))),
              reactions=("reaction", lambda s: " || ".join(uniq(s)))
          )
          .sort_values("gene"))

out.to_csv("atp_reactions_for_requested_genes.csv", index=False)
print(f"Wrote atp_reactions_for_requested_genes.csv with {len(out)} genes.")

# Diagnostic check for missing genes
all_input = {ln.strip().lower() for ln in open(GENES) if ln.strip()}
found = set(out['gene'].str.lower().unique())
missing = sorted(all_input - found)
pd.DataFrame(missing, columns=['missing_gene']).to_csv("missing_genes.csv", index=False)
print(len(missing), "genes not found. See missing_genes.csv for details.")

