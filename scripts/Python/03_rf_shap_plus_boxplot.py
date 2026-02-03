# -*- coding: utf-8 -*-
"""
Script 03 — Random Forest → SHAP on original protein features + OPTIONAL boxplot

Input:
- results/tables/gnn/GNN_input_matrix.xlsx

Outputs:
- results/figures/shap_protein/SHAP_global_bar_protein_RF.png
- results/figures/shap_protein/SHAP_heatmap_protein_RF.png
- results/figures/exploratory/boxplot_MWR_vs_MW_by_tag.png
Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
Author: Dr. El Hadji Malick Cisse
"""

import os
import random
import numpy as np
import pandas as pd

import shap
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier



# 0) CONFIG

IN_XLSX = os.path.join("results", "tables", "gnn", "GNN_input_matrix.xlsx")

OUT_SHAP = os.path.join("results", "figures", "shap_protein")
OUT_EXPL = os.path.join("results", "figures", "exploratory")

os.makedirs(OUT_SHAP, exist_ok=True)
os.makedirs(OUT_EXPL, exist_ok=True)

SEED = 42
MAKE_BOXPLOT = True  # <- set False if you don't want the exploratory figure

FEATURE_COLS = ["FC_vs_FW", "MC_vs_MW", "FWR_vs_FW", "MWR_vs_MW"]
TARGET_CLASSES = ["Waterlogging", "Recovery", "Waterlogging+SexDimorphism"]



# 1) Reproducibility

random.seed(SEED)
np.random.seed(SEED)



# 2) Tag simplification (same logic as Script 02)

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
    raise FileNotFoundError(f"Input not found: {IN_XLSX}")

df = pd.read_excel(IN_XLSX)

need = ["Tag"] + FEATURE_COLS
missing = [c for c in need if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in input: {missing}")

df["Simplified_Tag"] = df["Tag"].apply(simplify_tag)

# Focus on 3-class subset (as in your pipeline)
df_f = df[df["Simplified_Tag"].isin(TARGET_CLASSES)].copy()

X = df_f[FEATURE_COLS].astype(float)
y_labels = df_f["Simplified_Tag"].astype(str)

label_map = {lab: i for i, lab in enumerate(TARGET_CLASSES)}
y = y_labels.map(label_map).astype(int)



# 4) Random Forest classification

rf_model = RandomForestClassifier(
    n_estimators=200,
    max_depth=6,
    random_state=SEED
)
rf_model.fit(X, y)
print("Random Forest trained (protein features).")


# 5) SHAP on protein features (RF)

explainer = shap.TreeExplainer(rf_model)
shap_vals = explainer.shap_values(X)

# Normalize to list[class] -> array[n_samples, n_features]
if isinstance(shap_vals, np.ndarray):
    # sometimes (n_samples, n_features, n_classes)
    if shap_vals.ndim == 3 and shap_vals.shape[2] == len(TARGET_CLASSES):
        shap_list = [shap_vals[:, :, i] for i in range(len(TARGET_CLASSES))]
    elif shap_vals.ndim == 3 and shap_vals.shape[0] == len(TARGET_CLASSES):
        shap_list = [shap_vals[i, :, :] for i in range(len(TARGET_CLASSES))]
    else:
        raise ValueError(f"Unexpected SHAP array shape: {shap_vals.shape}")
else:
    shap_list = shap_vals

# 5.1 Global SHAP bar plot
plt.figure(figsize=(10, 4.5))
shap.summary_plot(
    shap_list,
    X,
    plot_type="bar",
    class_names=TARGET_CLASSES,
    show=False
)
plt.title("Global SHAP importance — protein features (Random Forest)", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(OUT_SHAP, "SHAP_global_bar_protein_RF.png"), dpi=300)
plt.close()

# 5.2 Heatmap across classes: mean(|SHAP|) per class × feature
heatmap_matrix = np.vstack([
    np.mean(np.abs(shap_list[i]), axis=0)
    for i in range(len(TARGET_CLASSES))
])

plt.figure(figsize=(8.5, 4.2))
sns.heatmap(
    heatmap_matrix,
    annot=True,
    fmt=".3f",
    cmap="YlGnBu",
    xticklabels=FEATURE_COLS,
    yticklabels=TARGET_CLASSES,
    cbar_kws={"label": "Mean(|SHAP value|)"}
)
plt.title("Protein feature importance by class (Random Forest + SHAP)")
plt.xlabel("Protein feature (contrast)")
plt.ylabel("Biological class")
plt.tight_layout()
plt.savefig(os.path.join(OUT_SHAP, "SHAP_heatmap_protein_RF.png"), dpi=300)
plt.close()

print("SHAP outputs saved in:", OUT_SHAP)


# 6) OPTIONAL: Exploratory boxplot (your original idea)

if MAKE_BOXPLOT:
    plt.figure(figsize=(8.5, 5))
    sns.boxplot(data=df, x="Simplified_Tag", y="MWR_vs_MW")
    plt.title("Distribution of MWR_vs_MW across biological tags")
    plt.xlabel("Biological tag")
    plt.ylabel("MWR_vs_MW (log2 spectrum merged)")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_EXPL, "boxplot_MWR_vs_MW_by_tag.png"), dpi=300)
    plt.close()
    print("Boxplot saved in:", OUT_EXPL)

print("Script 03 complete.")
