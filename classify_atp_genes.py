import json
import re
from collections import defaultdict

# Load the model file
with open("iML1515.json", "r") as f:
    M = json.load(f)

print("Model loaded:", M.get("id", "unknown"))
print("Number of reactions:", len(M.get("reactions", [])))
print("Number of genes:", len(M.get("genes", [])))

# Peek at a few reactions that involve ATP-like metabolites
print("\nChecking reactions that mention 'atp'...")
count = 0
for rxn in M["reactions"]:
    if any("atp" in m.lower() for m in rxn["metabolites"].keys()):
        print(rxn["id"], list(rxn["metabolites"].keys()))
        count += 1
        if count >= 5:
            break
# Identify reactions that consume ATP
rxn_to_atp = {}
rxn_genes = {}
for rxn in M["reactions"]:
    stoich = rxn.get("metabolites", {})
    consumes_atp = False
    if "atp_c" in stoich and "adp_c" in stoich:
        consumes_atp = True

    rxn_to_atp[rxn["id"]] = consumes_atp

    # Extract gene IDs from gene_reaction_rule
    gene_rule = rxn.get("gene_reaction_rule", "")
    genes = re.findall(r"b\d{4}", gene_rule)
    rxn_genes[rxn["id"]] = genes

# Count how many genes touch at least one ATP-consuming reaction
gene_atp = defaultdict(int)
for rxn_id, genes in rxn_genes.items():
    if rxn_to_atp[rxn_id]:
        for g in genes:
            gene_atp[g] += 1

print("\nGenes with ≥1 ATP-consuming reaction:",
      sum(1 for g in gene_atp if gene_atp[g] > 0))

# Save list of ATP-consuming genes to a text file
with open("atp_genes.txt", "w") as f:
    for g in sorted(gene_atp.keys()):
        if gene_atp[g] > 0:
            f.write(g + "\n")
print("Saved ATP-consuming gene list to atp_genes.txt")

# --- Map b-numbers to gene names using the model data ---
gene_id_to_name = {}
for g in M["genes"]:
    gene_id_to_name[g["id"]] = g.get("name", "")

# Create a CSV with both ID and gene name
with open("atp_genes_with_names.csv", "w") as f:
    f.write("b_number,gene_name\n")
    for g in sorted(gene_atp.keys()):
        if gene_atp[g] > 0:
            name = gene_id_to_name.get(g, "")
            f.write(f"{g},{name}\n")

print("✅ Saved ATP gene names to atp_genes_with_names.csv")

# --- Add reaction IDs for each gene ---
gene_to_rxns = defaultdict(list)
for rxn_id, genes in rxn_genes.items():
    if rxn_to_atp[rxn_id]:
        for g in genes:
            if gene_atp[g] > 0:
                gene_to_rxns[g].append(rxn_id)

# Save extended table with gene, name, and reactions
with open("atp_genes_with_reactions.csv", "w") as f:
    f.write("b_number,gene_name,reactions\n")
    for g in sorted(gene_atp.keys()):
        if gene_atp[g] > 0:
            name = gene_id_to_name.get(g, "")
            rxns = "|".join(sorted(set(gene_to_rxns[g])))
            f.write(f"{g},{name},{rxns}\n")

print("✅ Saved ATP gene → reaction mapping to atp_genes_with_reactions.csv")

# --- Add pathway (subsystem) info and save CSV ---

# Make a dictionary: reaction ID → subsystem name
rxn_to_pathway = {r["id"]: r.get("subsystem", "Unknown") for r in M["reactions"]}

# Make a dictionary: b-number → gene name
gene_to_name = {g["id"]: g.get("name", "") for g in M["genes"]}

# Write the output file
with open("atp_genes_with_pathways.csv", "w") as f:
    f.write("b_number,gene_name,reaction,subsystem\n")
    for rxn_id, genes in rxn_genes.items():
        if rxn_to_atp[rxn_id]:
            subsystem = rxn_to_pathway.get(rxn_id, "Unknown")
            for g in genes:
                f.write(f"{g},{gene_to_name.get(g, '')},{rxn_id},{subsystem}\n")

print("✅ Saved ATP gene–reaction–pathway list to atp_genes_with_pathways.csv")

# --- Condensed version: one line per gene ---
from collections import defaultdict

gene_summary = defaultdict(lambda: {"name": "", "reactions": set(), "subsystems": set()})

for rxn_id, genes in rxn_genes.items():
    if rxn_to_atp[rxn_id]:
        subsystem = rxn_to_pathway.get(rxn_id, "Unknown")
        for g in genes:
            gene_summary[g]["name"] = gene_to_name.get(g, "")
            gene_summary[g]["reactions"].add(rxn_id)
            gene_summary[g]["subsystems"].add(subsystem)

with open("atp_genes_condensed.csv", "w") as f:
    f.write("b_number,gene_name,reactions,subsystems\n")
    for g, d in gene_summary.items():
        f.write(f"{g},{d['name']},{'|'.join(sorted(d['reactions']))},{'|'.join(sorted(d['subsystems']))}\n")

print("✅ Saved condensed gene summary to atp_genes_condensed.csv")

