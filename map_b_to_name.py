import json, csv

with open("iML1515.json") as f:
    M = json.load(f)

with open("gene_mapping_name.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["b_number", "gene_name"])
    for g in M["genes"]:
        w.writerow([g["id"], g.get("name", "")])

print("✅ Saved b_number → gene_name map to gene_mapping_name.csv")

