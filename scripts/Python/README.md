Python scripts for proteomic sex dimorphic analysis with: matrix, classification, deep learning models, evaluation, and explainability.
# Build GNN-ready proteomics matrix
Run:
python scripts/Python/01_build_gnn_input_matrix.py

Inputs:
- data/gnn_inputs/*.xlsx (8 files: up/down for 4 contrasts; must contain Accession, Spectrum)
- data/taggedproteinsautoencoder.xlsx (optional; must contain Accession, Tag)

Outputs:
- results/tables/gnn/GNN_input_matrix.xlsx
# Script 02 — GCN → latent → XGBoost → SHAP
Run:
python scripts/Python/02_gcn_latent_xgb_shap.py

Input:
- results/tables/gnn/GNN_input_matrix.xlsx

Outputs:
- results/tables/gnn/GNN_latent_matrix.xlsx
- results/figures/shap_latent/*.png

# Script 03 — Random Forest → SHAP (protein features) + optional boxplot
Run:
python scripts/Python/03_rf_shap_plus_boxplot.py

Input:
- results/tables/gnn/GNN_input_matrix.xlsx

Outputs:
- results/figures/shap_protein/*.png
- (optional) results/figures/exploratory/boxplot_MWR_vs_MW_by_tag.png



