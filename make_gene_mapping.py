import json
import pandas as pd

with open("iML1515.json") as f:
    M = json.load(f)

rows = []
for g in M["genes"]:
    bnum = g["id"]
    locus = g.get("locus_tag", "")
    name = g.get("name", "")
    rows.append({"b_number": bnum, "locus_tag": locus, "gene_name": name})

pd.DataFrame(rows).to_csv("gene_mapping.csv", index=False)
print("✅ Saved mapping as gene_mapping.csv")

