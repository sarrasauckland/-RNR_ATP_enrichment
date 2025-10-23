import pandas as pd
from Bio import Entrez
import time

# --- Setup ---
Entrez.email = "your_email@example.com"  # replace with your email
batch_size = 10

# --- Load your SNP file ---
snps = pd.read_csv("mutated_genes.csv", sep="\t", dtype=str)
snps = snps[snps["locus_tag"].notna()]
unique_tags = snps["locus_tag"].unique().tolist()
print(f"Found {len(unique_tags)} unique locus tags to query...")

results = []

for i in range(0, len(unique_tags), batch_size):
    batch = unique_tags[i:i+batch_size]
    query = " OR ".join([f"{t}[Locus Tag]" for t in batch])
    print(f"Fetching {i+1}-{i+len(batch)} of {len(unique_tags)}...")

    try:
        handle = Entrez.esearch(db="gene", term=query, retmax=50)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            ids = record["IdList"]
            fetch_handle = Entrez.efetch(db="gene", id=",".join(ids), rettype="xml")
            records = Entrez.read(fetch_handle)
            fetch_handle.close()

            for r in records:
                locus_tag, bnum, name = None, None, None
                try:
                    name = r["Entrezgene_gene"]["Gene-ref"].get("Gene-ref_locus", None)
                    if "Gene-ref_locus-tag" in r["Entrezgene_gene"]["Gene-ref"]:
                        locus_tag = r["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus-tag"]
                    if "Entrezgene_comments" in r:
                        for comment in r["Entrezgene_comments"]:
                            if isinstance(comment, dict) and "Gene-commentary_text" in comment:
                                if comment["Gene-commentary_text"].startswith("b"):
                                    bnum = comment["Gene-commentary_text"]
                except Exception:
                    pass
                results.append({"locus_tag": locus_tag, "b_number": bnum, "gene_name": name})
        time.sleep(0.5)
    except Exception as e:
        print(f"⚠️ Error fetching batch {i}-{i+batch_size}: {e}")
        time.sleep(2)

# --- Save mapping ---
out = pd.DataFrame(results).drop_duplicates(subset=["locus_tag"])
out.to_csv("ecoli_locus_to_bnumber.csv", index=False)
print(f"✅ Saved {len(out)} mappings to ecoli_locus_to_bnumber.csv")

