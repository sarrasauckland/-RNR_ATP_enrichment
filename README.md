# -RNR_ATP_enrichment
Scripts and data for ATP-gene enrichment analyses in E. coli ΔRNR evolution experiments. Includes EcoCyc/BiGG ATP gene curation, SNP annotation, Fisher’s tests, and random-effects meta-analysis
| Script                              | Purpose                                                                                                          |
| :---------------------------------- | :--------------------------------------------------------------------------------------------------------------- |
| **`analyze_ecocyc_atp.py`**         | Parses EcoCyc reaction data to extract ATP-linked reactions (ATP, ADP, AMP, GTP substrates or products).         |
| **`filter_atp_reactions.py`**       | Filters out non-enzymatic, ambiguous, or indirect ATP reactions. Produces a curated list of ATP-consuming genes. |
| **`classify_atp_genes.py`**         | Cleans, merges, and classifies ATP genes by functional category (e.g., glycolysis, TCA, amino acid metabolism).  |
| **`compare_ecocyc_bigg.py`**        | Cross-references EcoCyc and BiGG models to validate and unify ATP gene annotations.                              |
| **`enrichment_ATP_fixed.py`**       | Performs Fisher’s exact tests for overall ATP-gene enrichment across all lines.                                  |
| **`linewise_enrichment_fixed.py`**  | Runs the same analysis per evolutionary line and exports summary tables.                                         |
| **`linewise_ATP_summary_fixed.py`** | Combines per-line outputs and performs DerSimonian–Laird meta-analysis for overall effect size.                  |
Python ≥ 3.9
Required libraries: pandas, numpy, scipy, statsmodels
