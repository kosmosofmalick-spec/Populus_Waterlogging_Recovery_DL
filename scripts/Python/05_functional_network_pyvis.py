# -*- coding: utf-8 -*-
"""
Script 05 — Interactive Functional Association Network (pyvis)

What it does
- Reads a curated table with Protein_ID, Fold_Change, COG and GO annotations
- Merges duplicated proteins by identical COG_Function_Description
- Builds a functional association network:
    - nodes = unique functions (COG_Function_Description)
    - edges = co-membership within the same COG functional category
- Exports:
    (1) Interactive HTML network (pyvis)
    (2) Cluster membership Excel
    (3) Cluster summary Excel

Input (expected columns)
- Protein_ID
- Fold_Change
- COG_ID
- COG_Function_Description
- COG_Functional_Categories
- GO_Biological_Process
- GO_Molecular_Function
- GO_Cellular_Component

Recommended repo layout
- results/tables/networks/Book.xlsx (file name)
- results/figures/networks/ (HTML output)
- results/tables/networks/ (cluster exports)
Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
Author: Dr. El Hadji Malick Cissse
"""

import os
from itertools import combinations
from collections import Counter

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
from pyvis.network import Network


# 0) CONFIG (edit only if needed)

INPUT_XLSX = os.path.join("results", "tables", "networks", "Book.xlsx")

OUT_DIR = os.path.join("results", "figures", "networks")
OUT_TABLES = os.path.join("results", "tables", "networks")

HTML_OUT = os.path.join(OUT_DIR, "PPI_FunctionalNetwork_Interactive.html")
CLUSTERS_OUT = os.path.join(OUT_TABLES, "PPI_FunctionalNetwork_Clusters.xlsx")
SUMMARY_OUT = os.path.join(OUT_TABLES, "PPI_Cluster_Summary.xlsx")

# Controls network density when categories are huge
MAX_EDGES_PER_CAT = 400
RANDOM_SEED = 42

# Visual scaling
DEFAULT_VMIN, DEFAULT_VMAX = 0.5, 2.5  #


# 1) IO checks

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(OUT_TABLES, exist_ok=True)

if not os.path.isfile(INPUT_XLSX):
    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_XLSX}\n"
        "Place your annotated table there, or edit INPUT_XLSX at the top."
    )


# 2) Load + validate

df = pd.read_excel(INPUT_XLSX)
df.columns = df.columns.astype(str).str.strip()

required = [
    "Protein_ID", "Fold_Change",
    "COG_ID", "COG_Function_Description", "COG_Functional_Categories",
    "GO_Biological_Process", "GO_Molecular_Function", "GO_Cellular_Component"
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Basic cleaning
df["COG_Function_Description"] = df["COG_Function_Description"].astype(str).str.strip()
df["COG_Functional_Categories"] = df["COG_Functional_Categories"].astype(str).str.strip()
df["Protein_ID"] = df["Protein_ID"].astype(str).str.strip()


# 3) Merge duplicates by FUNCTION DESCRIPTION

def combine_unique(series: pd.Series) -> str:
    vals = series.dropna().astype(str).str.strip().unique()
    vals = [v for v in vals if v != "" and v.lower() != "nan"]
    return "; ".join(vals) if len(vals) else np.nan

df_clean = (
    df.dropna(subset=["COG_Function_Description", "Fold_Change"])
      .groupby("COG_Function_Description", as_index=False)
      .agg({
          "Protein_ID": combine_unique,
          "Fold_Change": "mean",
          "COG_ID": combine_unique,
          "COG_Functional_Categories": combine_unique,
          "GO_Biological_Process": combine_unique,
          "GO_Molecular_Function": combine_unique,
          "GO_Cellular_Component": combine_unique
      })
)

print(f"Duplicates merged by function: {len(df)} → {len(df_clean)} unique functions")


# 4) Build edges based on shared COG functional categories

rng = np.random.default_rng(RANDOM_SEED)

edges = []
for cat, g in df_clean.groupby("COG_Functional_Categories", dropna=True):
    nodes = g["COG_Function_Description"].dropna().unique()
    if len(nodes) < 2:
        continue

    pairs = list(combinations(nodes, 2))

    # Downsample if too dense
    if len(pairs) > MAX_EDGES_PER_CAT:
        idx = rng.choice(len(pairs), size=MAX_EDGES_PER_CAT, replace=False)
        pairs = [pairs[i] for i in idx]

    edges.extend(pairs)

G = nx.Graph()
G.add_edges_from(edges)

# Ensure isolated nodes still appear
for fn in df_clean["COG_Function_Description"].dropna():
    if fn not in G:
        G.add_node(fn)

attrs = df_clean.set_index("COG_Function_Description").to_dict("index")
nx.set_node_attributes(G, attrs)

print(f"Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")


# 5) Visual encodings (Fold_Change → color & size)

fold_vals = np.array([G.nodes[n].get("Fold_Change", np.nan) for n in G.nodes()])
vmin, vmax = np.nanmin(fold_vals), np.nanmax(fold_vals)

if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or (vmin == vmax):
    vmin, vmax = DEFAULT_VMIN, DEFAULT_VMAX

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = mpl.colormaps.get_cmap("RdBu_r")  # blue ↔ white ↔ red

def fc_to_color(fc):
    if not np.isfinite(fc):
        fc = 1.0
    rgba = cmap(norm(fc))
    return mpl.colors.to_hex(rgba, keep_alpha=False)

def fc_to_size(fc):
    if not np.isfinite(fc):
        return 15.0
    # size grows with |fc - 1|
    denom = max(1e-9, abs(vmax - 1))
    return float(15 + 30 * abs(fc - 1) / denom)

def short_text(s, n=36):
    if not isinstance(s, str):
        return ""
    s = s.replace("\n", " ").strip()
    return s if len(s) <= n else s[: n - 3] + "..."


# 6) Create pyvis network

net = Network(
    height="900px",
    width="100%",
    bgcolor="#ffffff",
    font_color="#222222",
    directed=False,
    notebook=False
)

net.set_options("""
{
  "interaction": {"hover": true, "tooltipDelay": 120, "hideEdgesOnDrag": true},
  "nodes": {"borderWidth": 1, "shadow": {"enabled": true, "size": 3, "x":0, "y":1}},
  "edges": {"color": {"color":"#bbbbbb", "opacity":0.5}, "smooth": {"type": "continuous"}},
  "physics": {"stabilization": {"iterations": 200}}
}
""")


# 7) Add nodes + edges

for n in G.nodes():
    d = G.nodes[n]
    fc = d.get("Fold_Change", np.nan)

    label = short_text(n, 28)
    title = f"""
    <div style='font-family:Arial; font-size:12px;'>
      <b>Protein Function:</b> {n}<br>
      <b>Fold Change (mean):</b> {fc:.3f}<br>
      <b>Protein IDs:</b> {d.get('Protein_ID','')}<br>
      <b>COG ID:</b> {d.get('COG_ID','')}<br>
      <b>COG Category:</b> {d.get('COG_Functional_Categories','')}<br>
      <b>GO BP:</b> {d.get('GO_Biological_Process','')}<br>
      <b>GO MF:</b> {d.get('GO_Molecular_Function','')}<br>
      <b>GO CC:</b> {d.get('GO_Cellular_Component','')}<br>
    </div>
    """.replace("\n", "")

    net.add_node(
        n,
        label=label,
        title=title,
        color=fc_to_color(fc),
        size=fc_to_size(fc)
    )

for u, v in G.edges():
    net.add_edge(u, v, color="#cfcfcf", opacity=0.4)


# 8) Export interactive HTML

html_code = net.generate_html()

legend_html = """
<h3 style="font-family:Arial; color:#222; margin-top:8px; margin-bottom:2px;">
Protein Functional Association Network (merged by function)
</h3>
<p style="font-family:Arial; color:#444; margin-top:0; font-size:12px;">
Each node represents a unique protein function (COG_Function_Description). If multiple Protein IDs share
the same annotation, they are merged. Node color: fold change (blue→white→red). Hover nodes to view merged details.
Drag to rearrange; scroll to zoom.
</p>
"""

html_code = html_code.replace("<body>", f"<body>{legend_html}")

with open(HTML_OUT, "w", encoding="utf-8") as f:
    f.write(html_code)

print("Interactive HTML saved to:", HTML_OUT)


# 9) Cluster extraction + Excel exports

clusters = list(nx.connected_components(G))

cluster_data = []
for i, nodes in enumerate(clusters, start=1):
    sub = df_clean[df_clean["COG_Function_Description"].isin(nodes)].copy()
    sub["Cluster_ID"] = f"Cluster_{i}"
    cluster_data.append(sub)

clusters_df = pd.concat(cluster_data, ignore_index=True)

clusters_df = clusters_df[
    ["Cluster_ID", "COG_Function_Description", "Protein_ID", "Fold_Change",
     "COG_ID", "COG_Functional_Categories", "GO_Biological_Process",
     "GO_Molecular_Function", "GO_Cellular_Component"]
]

clusters_df.to_excel(CLUSTERS_OUT, index=False)
print("Cluster membership saved to:", CLUSTERS_OUT)

summary_rows = []
for i, nodes in enumerate(clusters, start=1):
    sub = df_clean[df_clean["COG_Function_Description"].isin(nodes)].copy()

    cog = sub["COG_Functional_Categories"].dropna().tolist()
    go_bp = sub["GO_Biological_Process"].dropna().tolist()

    main_cog = Counter(cog).most_common(1)[0][0] if cog else "-"
    main_bp = Counter(go_bp).most_common(1)[0][0] if go_bp else "-"

    summary_rows.append({
        "Cluster_ID": f"Cluster_{i}",
        "Function_Count": len(nodes),
        "Dominant_COG_Category": main_cog,
        "Dominant_GO_Biological_Process": main_bp
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_excel(SUMMARY_OUT, index=False)
print("Cluster summary saved to:", SUMMARY_OUT)
