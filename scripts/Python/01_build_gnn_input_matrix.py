# -*- coding: utf-8 -*-
"""
Build GNN-ready proteomics matrix from increased/decreased Excel files
- Reads 8 Excel files (4 contrasts × {increased, decreased})
- Keeps Accession + Spectrum, filters Spectrum > 0
- Computes log2(Spectrum)
- For each contrast, merges increased & decreased in proteins abundance and takes max(log2(Spectrum)) per protein
- Merges all contrasts into one wide matrix
- Optionally merges protein tags

Inputs:
- data/gnn_inputs/*.xlsx  (your increased/decreased in proteins abundance tables)
- data/taggedproteinsautoencoder.xlsx  (optional tagging table)

Outputs:
- results/tables/gnn/GNN_input_matrix.xlsx
Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
Author: Dr. El Hadji Malick Cisse
"""

import os
import pandas as pd
import numpy as np
from functools import reduce


# 0) CONFIG (edit these paths if your filenames differ)

BASE_DIR = os.path.join("data", "gnn_inputs")  # folder holding the 8 Excel files
TAG_TABLE = os.path.join("data", "taggedproteinsautoencoder.xlsx")  # optional
OUT_DIR = os.path.join("results", "tables", "gnn")
OUT_XLSX = os.path.join(OUT_DIR, "GNN_input_matrix.xlsx")

os.makedirs(OUT_DIR, exist_ok=True)

# Map of your 8 files (you can rename the keys, but keep structure)
file_info = {
    "FC12_vs_FW12_increased":   os.path.join(BASE_DIR, "FC12_vs_FW12_increased.xlsx"),
    "FC12_vs_FW12_decreased": os.path.join(BASE_DIR, "FC12_vs_FW12_decreased.xlsx"),
    "MC12_vs_MW12_increased":   os.path.join(BASE_DIR, "MC12_vs_MW12_increased.xlsx"),
    "MC12_vs_MW12_decreased": os.path.join(BASE_DIR, "MC12_vs_MW12_decreased.xlsx"),
    "FWR_vs_FW_increased":      os.path.join(BASE_DIR, "FWR_vs_FW_increased.xlsx"),
    "FWR_vs_FW_decreased":    os.path.join(BASE_DIR, "FWR_vs_FW_decreased.xlsx"),
    "MWR_vs_MW_increased":      os.path.join(BASE_DIR, "MWR_vs_MW_increased.xlsx"),
    "MWR_vs_MW_decreased":    os.path.join(BASE_DIR, "MWR_vs_MW_decreased.xlsx"),
}

# Group up/down into 4 contrasts
grouped_files = {
    "FC12_vs_FW12": ["FC12_vs_FW12_increased", "FC12_vs_FW12_decreased"],
    "MC12_vs_MW12": ["MC12_vs_MW12_increased", "MC12_vs_MW12_decreased"],
    "FWR_vs_FW":    ["FWR_vs_FW_increased", "FWR_vs_FW_decreased"],
    "MWR_vs_MW":    ["MWR_vs_MW_increased", "MWR_vs_MW_decreased"],
}


# 1) Helper: read one file safely and compute log2 spectrum

def read_and_prepare(path: str, out_col: str) -> pd.DataFrame:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing input file: {path}")

    df = pd.read_excel(path)

    needed = {"Accession", "Spectrum"}
    if not needed.issubset(set(df.columns)):
        raise ValueError(
            f"File {os.path.basename(path)} must contain columns {needed}. "
            f"Found: {list(df.columns)}"
        )

    df = df[["Accession", "Spectrum"]].dropna(subset=["Spectrum"])
    df = df[df["Spectrum"] > 0].copy()

    df[out_col] = np.log2(df["Spectrum"].astype(float))
    df = df[["Accession", out_col]]

    # If Accession duplicates exist, keep the max log2_Spectrum
    df = df.groupby("Accession", as_index=False)[out_col].max()

    return df



# 2) Build per-contrast matrix = max(increased, decreased)

merged_data = {}

for contrast, parts in grouped_files.items():
    up_key, down_key = parts[0], parts[1]

    df_up = read_and_prepare(file_info[up_key], out_col=f"{contrast}_increased")
    df_dn = read_and_prepare(file_info[down_key], out_col=f"{contrast}_decreased")

    temp = pd.merge(df_up, df_dn, on="Accession", how="outer")

    # Take max of increased/dcreeased; missing becomes NaN then filled later
    temp[contrast] = temp[[f"{contrast}_increased", f"{contrast}_decreased"]].max(axis=1)

    merged_data[contrast] = temp[["Accession", contrast]]


# 3) Merge all contrasts into final matrix

final_merged = reduce(
    lambda left, right: pd.merge(left, right, on="Accession", how="outer"),
    merged_data.values(),
)

final_merged = final_merged.fillna(0)


# 4) Optional: merge protein tags

if os.path.isfile(TAG_TABLE):
    tagging_df = pd.read_excel(TAG_TABLE)

    if "Accession" not in tagging_df.columns:
        raise ValueError("Tag table must contain column: Accession")

    # If multiple tag columns exist, keep all; if Tag missing, set None
    final_merged = pd.merge(final_merged, tagging_df, on="Accession", how="left")
    if "Tag" in final_merged.columns:
        final_merged["Tag"] = final_merged["Tag"].fillna("None")
else:
    print(f"NOTE: Tag table not found at {TAG_TABLE}. Skipping Tag merge.")

# 5) Save output

final_merged.to_excel(OUT_XLSX, index=False)
print(" GNN input matrix created successfully!")
print("Saved:", OUT_XLSX)
print("Shape:", final_merged.shape)
