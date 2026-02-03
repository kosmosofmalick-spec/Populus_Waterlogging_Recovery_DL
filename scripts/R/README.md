R scripts for proteomics sex dimorphic analysis, normalization, differential abundance analysis, enrichment, and plotting.
# UpSet + Venn (DEG set overlaps)
Run:
Rscript scripts/R/01_upset_venn_multiset.R

Inputs:
- data/Increased_proteins_abundance.xlsx
- data/Decreased_proteins_abundance.xlsx

Outputs:
- results/figures/upset_venn/*.png
# GO enrichment visualization (proteins)
Run:
Rscript scripts/R/02_GO_enrichment_visualization_proteins.R

Input:
- data/GO_enrichment.xlsx

Outputs:
- results/figures/go_enrichment/*.pdf and *.png
# Latent screening + sex-dimorphism heatmaps (proteins)
Run:
Rscript scripts/R/03_latent_screening_sex_dimorphism_proteins.R

Input:
- data/Proteomics_with_Latent.xlsx

Outputs:
- results/tables/latent_screening/Selected_DECREASED_*.xlsx
- results/tables/latent_screening/Selected_INCREASED_*.xlsx
- results/figures/heatmaps/latent_screening/Heatmap_*.png

# Sex dimorphism figures (morpho-physiological-scale-free)
Run:
Rscript scripts/R/04_sex_dimorphism_scale_free_figures.R

Requires outputs from:
- scripts/Python/04_sex_dimorphism_scale_free_contrasts.py

Outputs:
- results/figures/sex_dimorphism/*.png







