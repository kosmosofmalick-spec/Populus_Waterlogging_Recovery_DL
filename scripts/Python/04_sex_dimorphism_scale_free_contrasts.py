# -*- coding: utf-8 -*-
"""
Morpho-physiological Sex dimorphism & recovery contrasts from Treatment×Sex means (scale-free)

Inputs:
- data/Data_to_analyse.xlsx

Expects columns (case-sensitive in file):
- 'Treatments', 'Sex' + trait columns (e.g., Shoot height, Basal Stem diameter, Pn, Wue, ...)

Notes:
- This script strips whitespace from column names (e.g., 'Pn ' -> 'Pn').

Outputs:
- results/tables/sex_dimorphism/sex_dimorphism_outputs.xlsx (multi-sheet)
- results/tables/sex_dimorphism/*.csv (all contrast tables)
- results/tables/sex_dimorphism/traits_means_clean.csv (cleaned copy with stripped headers)
Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
Author: Dr. El Hadji Malick Cisse
"""

import os
import sys
import pandas as pd
import numpy as np


# 0) CONFIG

INPUT_XLSX = os.path.join("data", "Data_to_analyse.xlsx")

OUT_DIR = os.path.join("results", "tables", "sex_dimorphism")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_EXCEL = os.path.join(OUT_DIR, "sex_dimorphism_outputs.xlsx")


# 1) Helpers

def zscore_cols(df, cols):
    """Z-score standardization across all rows for the given columns."""
    Z = df.copy()
    for c in cols:
        mu = Z[c].mean()
        sd = Z[c].std(ddof=0)
        Z[c] = 0.0 if (sd == 0 or pd.isna(sd)) else (Z[c] - mu) / sd
    return Z

def compute_sex_diffs_by_treatment(pivot, treatments, trait_cols, sexA="Male", sexB="Female"):
    """For each treatment, compute (sexA - sexB) per trait."""
    out = {}
    for t in treatments:
        row = {tr: pivot.loc[t, (tr, sexA)] - pivot.loc[t, (tr, sexB)] for tr in trait_cols}
        out[t] = row
    df = pd.DataFrame.from_dict(out, orient="index")
    df.index.name = "Treatment"
    return df

def compute_within_sex_contrast(pivot, trait_cols, sex, treatA, treatB):
    """Within-sex contrast (treatA - treatB) across traits."""
    row = {tr: pivot.loc[treatA, (tr, sex)] - pivot.loc[treatB, (tr, sex)] for tr in trait_cols}
    return pd.Series(row, name=f"{sex}:{treatA}-vs-{treatB}")

def contribution_table_from_diffs(diffs_like):
    """
    Compute mean absolute standardized difference per trait and normalize to % contribution.
    """
    if isinstance(diffs_like, pd.Series):
        mean_abs = diffs_like.abs()
    else:
        mean_abs = diffs_like.abs().mean(axis=0)

    total = mean_abs.sum()
    pct = mean_abs * 0.0 if (total == 0 or pd.isna(total)) else (100.0 * mean_abs / total)

    return pd.DataFrame(
        {"mean_abs_std_diff": mean_abs.values, "percent_contribution": pct.values},
        index=mean_abs.index
    )

def save_series_or_df(obj, out_path):
    if isinstance(obj, pd.Series):
        obj.to_frame().T.to_csv(out_path, index=True)
    else:
        obj.to_csv(out_path, index=True)

# 2) Load & validate

if not os.path.isfile(INPUT_XLSX):
    sys.exit(f"ERROR: Input file not found:\n{INPUT_XLSX}")

df_raw = pd.read_excel(INPUT_XLSX, sheet_name=0)

# strip whitespace from column names
df = df_raw.copy()
df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]

REQ_ID = ["Treatments", "Sex"]
missing_ids = [c for c in REQ_ID if c not in df.columns]
if missing_ids:
    sys.exit(f"ERROR: Required columns missing: {missing_ids}\nFound: {df.columns.tolist()}")

trait_cols = [c for c in df.columns if c not in REQ_ID]
if not trait_cols:
    sys.exit("ERROR: No trait columns found (besides 'Treatments' and 'Sex').")

treatments = df["Treatments"].astype(str).unique().tolist()
sexes = df["Sex"].astype(str).unique().tolist()

# 3) Standardize and build pivot

Z = zscore_cols(df, trait_cols)
pivot = Z.pivot(index="Treatments", columns="Sex", values=trait_cols).sort_index()

# 4) Contrasts

sex_diffs = compute_sex_diffs_by_treatment(pivot, treatments, trait_cols, sexA="Male", sexB="Female")

have_CR = {"Control", "Recovery"}.issubset(set(treatments))
cr_m = compute_within_sex_contrast(pivot, trait_cols, "Male", "Control", "Recovery") if have_CR else None
cr_f = compute_within_sex_contrast(pivot, trait_cols, "Female", "Control", "Recovery") if have_CR else None

have_WR = {"Waterlogging", "Recovery"}.issubset(set(treatments))
wr_m = compute_within_sex_contrast(pivot, trait_cols, "Male", "Waterlogging", "Recovery") if have_WR else None
wr_f = compute_within_sex_contrast(pivot, trait_cols, "Female", "Waterlogging", "Recovery") if have_WR else None

cr_inter = (cr_m - cr_f).rename("CrossSex_Interaction:Control-Rec") if (cr_m is not None and cr_f is not None) else None
wr_inter = (wr_m - wr_f).rename("CrossSex_Interaction:Waterlog-Rec") if (wr_m is not None and wr_f is not None) else None

# 5) Percent contributions

contrib_frames = []
contrib_sex = contribution_table_from_diffs(sex_diffs); contrib_sex["contrast"] = "Sex(M−F)|byTreatment"; contrib_frames.append(contrib_sex)

if cr_m is not None:
    cm = contribution_table_from_diffs(cr_m); cm["contrast"] = "Male:Control−Recovery"; contrib_frames.append(cm)
if cr_f is not None:
    cf = contribution_table_from_diffs(cr_f); cf["contrast"] = "Female:Control−Recovery"; contrib_frames.append(cf)
if wr_m is not None:
    wm = contribution_table_from_diffs(wr_m); wm["contrast"] = "Male:Waterlogging−Recovery"; contrib_frames.append(wm)
if wr_f is not None:
    wf = contribution_table_from_diffs(wr_f); wf["contrast"] = "Female:Waterlogging−Recovery"; contrib_frames.append(wf)
if cr_inter is not None:
    ci = contribution_table_from_diffs(cr_inter); ci["contrast"] = "CrossSex:(C−R)"; contrib_frames.append(ci)
if wr_inter is not None:
    wi = contribution_table_from_diffs(wr_inter); wi["contrast"] = "CrossSex:(W−R)"; contrib_frames.append(wi)

contrib_all = pd.concat(contrib_frames, axis=0)
contrib_all.index.name = "Trait"
contrib_all = contrib_all.reset_index()[["Trait", "mean_abs_std_diff", "percent_contribution", "contrast"]]

# 6) Save outputs

csv_paths = {}

def save(obj, name):
    path = os.path.join(OUT_DIR, name)
    save_series_or_df(obj, path)
    csv_paths[name] = path

save(sex_diffs, "Z_sex_diffs_by_treatment.csv")
if cr_m is not None: save(cr_m, "Z_cr_m.csv")
if cr_f is not None: save(cr_f, "Z_cr_f.csv")
if wr_m is not None: save(wr_m, "Z_wr_m.csv")
if wr_f is not None: save(wr_f, "Z_wr_f.csv")
if cr_inter is not None: save(cr_inter, "Z_cr_inter.csv")
if wr_inter is not None: save(wr_inter, "Z_wr_inter.csv")

contrib_csv = os.path.join(OUT_DIR, "Z_contributions_all_contrasts.csv")
contrib_all.to_csv(contrib_csv, index=False)
csv_paths["Z_contributions_all_contrasts.csv"] = contrib_csv

clean_means_csv = os.path.join(OUT_DIR, "traits_means_clean.csv")
df.to_csv(clean_means_csv, index=False)
csv_paths["traits_means_clean.csv"] = clean_means_csv

with pd.ExcelWriter(OUT_EXCEL, engine="xlsxwriter") as xw:
    df_raw.to_excel(xw, sheet_name="Original_Data", index=False)
    df.to_excel(xw, sheet_name="Cleaned_Data", index=False)
    Z.to_excel(xw, sheet_name="Z_Scored_Data", index=False)

    sex_diffs.to_excel(xw, sheet_name="Z_SexDiffs_byTreatment")
    if cr_m is not None: cr_m.to_frame().T.to_excel(xw, sheet_name="Z_CR_Male")
    if cr_f is not None: cr_f.to_frame().T.to_excel(xw, sheet_name="Z_CR_Female")
    if wr_m is not None: wr_m.to_frame().T.to_excel(xw, sheet_name="Z_WR_Male")
    if wr_f is not None: wr_f.to_frame().T.to_excel(xw, sheet_name="Z_WR_Female")
    if cr_inter is not None: cr_inter.to_frame().T.to_excel(xw, sheet_name="Z_Interaction_CR")
    if wr_inter is not None: wr_inter.to_frame().T.to_excel(xw, sheet_name="Z_Interaction_WR")

    contrib_all.to_excel(xw, sheet_name="Z_Contrib_AllContrasts", index=False)

print("\nDone. Outputs written to:", OUT_DIR)
print("Main workbook:", OUT_EXCEL)
print("CSVs:")
for name, path in csv_paths.items():
    print(" -", path)
