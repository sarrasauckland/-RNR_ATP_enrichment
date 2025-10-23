import pandas as pd
import re

def _read_csv_flexible(path):
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def _split_multi(val):
    if not isinstance(val, str) or val.strip() == '':
        return []
    # Split on commas, semicolons, pipes, or single/double slashes
    parts = re.split(r'[;,|/]{1,2}\s*', val.strip())
    return [p.strip() for p in parts if p.strip()]

def analyze(file_path):
    df = _read_csv_flexible(file_path)
    col_substrates = next((c for c in df.columns if "substrate" in c), None)
    col_reaction = next((c for c in df.columns if "reaction" in c), None)
    col_enzymes = next((c for c in df.columns if "enzyme" in c), None)
    col_genes = next((c for c in df.columns if "gene" in c), None)

    if not col_substrates:
        raise ValueError("No substrates column found.")

    filtered = df[df[col_substrates].str.upper().str.contains("ATP|GTP|NTP", na=False)]

    records = []
    for _, row in filtered.iterrows():
        rxn = row.get(col_reaction, "unknown")
        enzymes = _split_multi(row.get(col_enzymes, ""))
        genes = _split_multi(row.get(col_genes, ""))
        for g in genes or [""]:
            for e in enzymes or [""]:
                records.append(dict(reaction=rxn, enzyme=e, gene=g))

    out = pd.DataFrame(records)
    out.to_csv("reaction_gene_map_atp.tsv", sep="\t", index=False)
    gene_counts = out[out['gene'] != ""].groupby('gene')['reaction'].nunique().reset_index()
    gene_counts.columns = ['gene', 'n_atp_consuming_reactions']
    gene_counts.to_csv("gene_counts_atp.tsv", sep="\t", index=False)
    print(f"Wrote {len(out)} reaction–gene pairs and {len(gene_counts)} unique genes.")

if __name__ == "__main__":
    analyze("ecocyc_all_reactions.csv")

