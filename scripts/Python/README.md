Python scripts for proteomic sex dimorphic analysis with: matrix, classification, deep learning models, evaluation, and explainability.
# Build GNN-ready proteomics matrix
Run:
python scripts/Python/01_build_gnn_input_matrix.py

Inputs:
- data/gnn_inputs/*.xlsx (8 files: up/down for 4 contrasts; must contain Accession, Spectrum)
- data/taggedproteinsautoencoder.xlsx (optional; must contain Accession, Tag)

Outputs:
- results/tables/gnn/GNN_input_matrix.xlsx


