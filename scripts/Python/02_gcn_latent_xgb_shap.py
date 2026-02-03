# -*- coding: utf-8 -*-
"""
Script 02 — GCN → latent embedding → XGBoost → SHAP
Populus deltoides proteomics (waterlogging/recovery/sex dimorphism)

Input:
- results/tables/gnn/GNN_input_matrix.xlsx  (from Script 01)

Outputs:
- results/tables/gnn/GNN_latent_matrix.xlsx
- results/figures/shap_latent/SHAP_global_bar_latent.png
- results/figures/shap_latent/SHAP_beeswarm_latent_<class>.png
- results/figures/shap_latent/SHAP_top_features_latent_<class>.png

Notes:
- Graph topology: KNN in feature space (k=5). Data-driven; no external PPI required.
- Reproducible seeds are fixed.
Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
Author: Dr. El Hadji Malick Cisse
"""

import os
import random
import numpy as np
import pandas as pd

import torch
import torch.nn.functional as F

import xgboost as xgb
import shap

import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder
from sklearn.neighbors import NearestNeighbors

from torch_geometric.data import Data
from torch_geometric.nn import GCNConv



# 0) CONFIG

IN_XLSX = os.path.join("results", "tables", "gnn", "GNN_input_matrix.xlsx")

OUT_TABLES = os.path.join("results", "tables", "gnn")
OUT_FIGS = os.path.join("results", "figures", "shap_latent")

OUT_LATENT_XLSX = os.path.join(OUT_TABLES, "GNN_latent_matrix.xlsx")

os.makedirs(OUT_TABLES, exist_ok=True)
os.makedirs(OUT_FIGS, exist_ok=True)

SEED = 42
K_NEIGHBORS = 5
EPOCHS = 150
LR = 0.01
WEIGHT_DECAY = 5e-4

FEATURE_COLS = ["FC_vs_FW", "MC_vs_MW", "FWR_vs_FW", "MWR_vs_MW"]
TARGET_CLASSES = ["Waterlogging", "Recovery", "Waterlogging+SexDimorphism"]



# 1) Reproducibility

def set_seed(seed: int = 42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

set_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False



# 2) Tag simplification

def simplify_tag(tag):
    tag = str(tag)
    if ("Waterlogging" in tag) and ("Sex_Dimorphism" in tag):
        return "Waterlogging+SexDimorphism"
    elif "Waterlogging" in tag:
        return "Waterlogging"
    elif "Recovery" in tag:
        return "Recovery"
    elif "Sex_Dimorphism" in tag:
        return "Sex_Dimorphism"
    else:
        return "None"



# 3) Load data

if not os.path.isfile(IN_XLSX):
    raise FileNotFoundError(
        f"Input not found: {IN_XLSX}\n"
        "Run scripts/Python/01_build_gnn_input_matrix.py first."
    )

df = pd.read_excel(IN_XLSX)

need = ["Accession", "Tag"] + FEATURE_COLS
missing = [c for c in need if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in input: {missing}")

df["Simplified_Tag"] = df["Tag"].apply(simplify_tag)

X = df[FEATURE_COLS].to_numpy(dtype=float)
accessions = df["Accession"].astype(str).tolist()

label_encoder = LabelEncoder()
y_all = label_encoder.fit_transform(df["Simplified_Tag"].astype(str))



# 4) Build KNN graph (k=5) in feature space

nbrs = NearestNeighbors(n_neighbors=min(K_NEIGHBORS + 1, len(df)), metric="euclidean")
nbrs.fit(X)
_, idxs = nbrs.kneighbors(X)

edges = set()
for i in range(idxs.shape[0]):
    for j in idxs[i]:
        if i == j:
            continue
        a, b = (i, j) if i < j else (j, i)
        edges.add((a, b))

edge_index = torch.tensor(list(edges), dtype=torch.long).t().contiguous()
if edge_index.numel() == 0:
    raise RuntimeError("Graph has no edges. Check K_NEIGHBORS or input size.")



# 5) GCN model

x_tensor = torch.tensor(X, dtype=torch.float)
y_tensor = torch.tensor(y_all, dtype=torch.long)
data = Data(x=x_tensor, edge_index=edge_index, y=y_tensor)

class GCN(torch.nn.Module):
    def __init__(self, in_dim: int, hidden: int = 8, out_dim: int = 4):
        super().__init__()
        self.conv1 = GCNConv(in_dim, hidden)
        self.conv2 = GCNConv(hidden, out_dim)

    def forward(self, d):
        x, edge_index = d.x, d.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training, p=0.3)
        x = self.conv2(x, edge_index)
        return x

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = GCN(in_dim=len(FEATURE_COLS), hidden=8, out_dim=4).to(device)
data = data.to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)



# 6) Train GCN (supervised on Simplified_Tag)

model.train()
for epoch in range(EPOCHS):
    optimizer.zero_grad()
    out = model(data)
    loss = F.cross_entropy(out, data.y)
    loss.backward()
    optimizer.step()

    if epoch % 30 == 0:
        print(f"Epoch {epoch:03d} | Loss: {loss.item():.4f}")



# 7) Extract latent features (3 dims)

model.eval()
with torch.no_grad():
    latent = model(data).detach().cpu().numpy()

df["Latent_1"] = latent[:, 0]
df["Latent_2"] = latent[:, 1]
df["Latent_3"] = latent[:, 2]

df.to_excel(OUT_LATENT_XLSX, index=False)
print("✅ Latent matrix saved:", OUT_LATENT_XLSX)


# ---------------------------------------------------------
# 8) XGBoost on latent space (2D for interpretability)
# ---------------------------------------------------------
df_f = df[df["Simplified_Tag"].isin(TARGET_CLASSES)].copy()
X_latent = df_f[["Latent_1", "Latent_2"]]
y_labels = df_f["Simplified_Tag"].astype(str)

label_map = {lab: i for i, lab in enumerate(TARGET_CLASSES)}
y_enc = y_labels.map(label_map).astype(int)

xgb_model = xgb.XGBClassifier(
    objective="multi:softprob",
    num_class=len(TARGET_CLASSES),
    n_estimators=100,
    max_depth=4,
    learning_rate=0.1,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=SEED,
    eval_metric="mlogloss",
)
xgb_model.fit(X_latent, y_enc)
print("✅ XGBoost trained on latent space.")


# 9) SHAP on latent space

explainer = shap.TreeExplainer(xgb_model)
shap_vals = explainer.shap_values(X_latent)

# Normalize SHAP output to list[class] -> array[n_samples, n_features]
if isinstance(shap_vals, np.ndarray):
    if shap_vals.ndim == 3 and shap_vals.shape[2] == len(TARGET_CLASSES):
        shap_list = [shap_vals[:, :, i] for i in range(len(TARGET_CLASSES))]
    elif shap_vals.ndim == 3 and shap_vals.shape[0] == len(TARGET_CLASSES):
        shap_list = [shap_vals[i, :, :] for i in range(len(TARGET_CLASSES))]
    else:
        raise ValueError(f"Unexpected SHAP array shape: {shap_vals.shape}")
else:
    shap_list = shap_vals

# 9.1 Global bar plot
plt.figure(figsize=(10, 4.5))
shap.summary_plot(
    shap_list,
    X_latent,
    plot_type="bar",
    class_names=TARGET_CLASSES,
    show=False
)
plt.title("Global SHAP importance — latent variables (XGBoost)", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUT_FIGS, "SHAP_global_bar_latent.png"), dpi=300)
plt.close()

# 9.2 Per-class bar + beeswarm
for i, lab in enumerate(TARGET_CLASSES):
    vals = np.array(shap_list[i])

    # Per-class mean |SHAP| bar (simple Matplotlib)
    mean_abs = np.abs(vals).mean(axis=0)
    order = np.argsort(mean_abs)[::-1]
    feat_names = list(X_latent.columns)

    plt.figure(figsize=(6.5, 4.2))
    plt.barh([feat_names[j] for j in order], mean_abs[order])
    plt.gca().invert_yaxis()
    plt.xlabel("Mean |SHAP value|")
    plt.title(f"Top SHAP features (latent) — {lab}")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_FIGS, f"SHAP_top_features_latent_{lab}.png"), dpi=300)
    plt.close()

    # Beeswarm
    plt.figure(figsize=(8, 5))
    shap.summary_plot(vals, X_latent, plot_type="dot", show=False)
    plt.title(f"SHAP beeswarm (latent) — {lab}", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_FIGS, f"SHAP_beeswarm_latent_{lab}.png"), dpi=300)
    plt.close()

print("Script 02 complete. Latent + SHAP figures saved in:", OUT_FIGS)
