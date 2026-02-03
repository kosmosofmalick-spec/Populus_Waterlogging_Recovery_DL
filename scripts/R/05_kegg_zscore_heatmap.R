# ==========================================================
# KEGG Pathway Z-score Heatmap (Excel version, styled)
# Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
# Input: long table with columns: FileName, Variable, Zscore
# Output: PDF + PNG heatmap
# Author: Dr. El Hadji Malick Cisse
# ==========================================================

  library(tidyverse)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(grid)

# 1) Paths

input_path <- file.path("data", "KEGG_Zscore.xlsx") 
out_dir <- file.path("results", "figures", "kegg")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop(paste0(
    "Input file not found: ", input_path, "\n",
    "Please place your Excel file in /data or edit `input_path`."
  ))
}

# 2) Read Excel data

# Sheet 1 by default
df <- readxl::read_excel(input_path, sheet = 1)

# 3) Validate required columns

req <- c("FileName", "Variable", "Zscore")
if (!all(req %in% names(df))) {
  stop("Excel file must contain columns: FileName, Variable, and Zscore.")
}

# 4) Summarize duplicates (mean Zscore per FileName × Variable)

df_sum <- df %>%
  group_by(FileName, Variable) %>%
  summarise(Zscore = mean(Zscore, na.rm = TRUE), .groups = "drop")

# 5) Transform to wide matrix

mat_full <- df_sum %>%
  pivot_wider(names_from = FileName, values_from = Zscore) %>%
  column_to_rownames("Variable") %>%
  as.matrix()

# 6) Keep only pathways present in all comparisons

mat_inter <- mat_full[complete.cases(mat_full), , drop = FALSE]
if (nrow(mat_inter) == 0) stop("No complete-case pathways found across all comparisons.")

# 7) Select top 13 pathways by variability (row SD)

row_sd <- matrixStats::rowSds(mat_inter, na.rm = TRUE)
ord <- order(row_sd, decreasing = TRUE)
mat_top <- mat_inter[ord[seq_len(min(13, nrow(mat_inter)))], , drop = FALSE]

# 8) Color scaling

lim <- max(abs(mat_top), na.rm = TRUE)
col_fun <- circlize::colorRamp2(
  seq(-lim, lim, length.out = 7),
  c("#00441b", "#1b7837", "#a6d96a", "white",
    "#fdae61", "#d7191c", "#67001f")
)

# 9) Build heatmap

ht <- Heatmap(
  mat_top,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Z-score",
    legend_direction = "horizontal",
    legend_width = unit(5, "cm"),
    title_gp = gpar(fontface = "bold")
  ),
  column_title = "Comparisons",
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  row_title = "KEGG pathway",
  row_title_gp = gpar(fontface = "bold", fontsize = 12)
)

# 10) Save PDF + PNG

pdf_out <- file.path(out_dir, "KEGG_Zscore_Heatmap.pdf")
png_out <- file.path(out_dir, "KEGG_Zscore_Heatmap.png")

pdf(pdf_out, width = 7, height = 8)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

png(png_out, width = 2200, height = 2600, res = 300)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

message(" Saved: ", pdf_out)
message(" Saved: ", png_out)
