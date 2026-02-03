# Populus_Waterlogging_Recovery_DL
R and Python scripts for deep learning-assisted proteomics analysis of Populus under waterlogging and recovery
# Deep learning-assisted proteomic analysis of Populus deltoides under waterlogging and recovery

This repository contains the analysis scripts used in the manuscript:
**[Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and subsequent recovery]**

## Contents
- `scripts/R/` — R scripts for enrichment plots, heatmaps, UpSet/Venn, latent screening visualization, KEGG Z-score heatmap
- `scripts/Python/` — Python scripts for matrix assembly, GCN→latent→XGBoost→SHAP, RF→SHAP (+ optional boxplots), and interactive functional association network (pyvis)
- `scripts/` READMEs — how to run each script and required inputs/outputs

## Data
All raw and processed datasets used by these scripts are deposited on Zenodo:
- Zenodo record: **[(https://doi.org/10.5281/zenodo.17491527)]**

## How to run
See `scripts/README.md` and the READMEs inside `scripts/R/` and `scripts/Python/`.

## License
MIT for code
