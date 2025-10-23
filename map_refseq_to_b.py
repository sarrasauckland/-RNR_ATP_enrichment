import json
import csv

# Load iML1515 model
with open("iML1515.json") as f:
    M = json.load(f)

# Write out a mapping file of b_number <-> locus_tag
with open("gene_mapping.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["b_number", "locus_tag"])
    for g in M["genes"]:
        w.writerow([g["id"], g.get("locus_tag", "")])

print("✅ Saved mapping to gene_mapping.csv")

