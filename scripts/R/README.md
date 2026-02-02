R scripts for proteomics sex dimorphic analysis, normalization, differential abundance analysis, enrichment, and plotting.
# UpSet + Venn (DEG set overlaps)
Run:
Rscript scripts/R/01_upset_venn_multiset.R

Inputs:
- data/Increased_proteins_abundance.xlsx
- data/Decreased_proteins_abundance.xlsx

Outputs:
- results/figures/upset_venn/*.png




