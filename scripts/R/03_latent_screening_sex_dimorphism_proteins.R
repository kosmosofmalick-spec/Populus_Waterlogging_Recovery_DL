# ==========================================================
# Sex-dimorphism latent screening & visualization (proteomics)
# Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
# - Select top DECREASED and INCREASED proteins by latent magnitude
# - Export selected tables
# - Generate heatmaps with latent magnitude annotations
# Author: Dr. El Hadji Malick Cisse
# ==========================================================

  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)


# 1) Paths

in_file <- file.path("data", "Proteomics_with_Latent.xlsx")

out_tables <- file.path("results", "tables", "latent_screening")
out_figs   <- file.path("results", "figures", "heatmaps", "latent_screening")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop(paste0(
    "Input file not found: ", in_file, "\n",
    "This repo is code-only. Place Proteomics_with_Latent.xlsx in /data\n",
    "or edit `in_file` to your local path."
  ))
}

today <- format(Sys.Date(), "%Y-%m-%d")

# 2) Load data

df <- readxl::read_excel(in_file)

# Required columns check
req_cols <- c("Tag", "Latent_1", "Latent_2", "Latent_3", "Accession")
missing <- setdiff(req_cols, names(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# 3) Keep only sex-dimorphism proteins (waterlogging ± recovery)

df_sd <- df %>%
  filter(Tag %in% c("Waterlogging;Sex_Dimorphism",
                    "Waterlogging;Sex_Dimorphism;Recovery"))

# 4) Latent magnitude

df_sd <- df_sd %>%
  mutate(Latent_Magnitude = sqrt(Latent_1^2 + Latent_2^2 + Latent_3^2))

# 5) Identify FC columns & split by direction labels in column names

fc_cols <- grep("^FC_", names(df_sd), value = TRUE)

fc_decreased_cols <- grep("down", fc_cols, value = TRUE, ignore.case = TRUE)
fc_increased_cols <- grep("up",   fc_cols, value = TRUE, ignore.case = TRUE)

if (length(fc_decreased_cols) == 0 && length(fc_increased_cols) == 0) {
  stop("No FC columns detected containing 'up' or 'down'. Please check your column names.")
}

# 6) Classification: DECREASED vs INCREASED vs MIXED

classify_by_fc <- function(row) {
  vals_dec <- as.numeric(unlist(row[fc_decreased_cols]))
  vals_inc <- as.numeric(unlist(row[fc_increased_cols]))

  has_dec <- any(is.finite(vals_dec) & !is.na(vals_dec))
  has_inc <- any(is.finite(vals_inc) & !is.na(vals_inc))

  if (has_dec && !has_inc) return("DECREASED")
  if (has_inc && !has_dec) return("INCREASED")
  return("MIXED")
}

df_sd$Direction <- apply(df_sd, 1, classify_by_fc)

# 7) Filter proteins with enough valid FC values per group

min_valid <- 2
n_take <- 30

decreased_set <- df_sd %>%
  filter(Direction == "DECREASED") %>%
  filter(rowSums(!is.na(select(., all_of(fc_decreased_cols)))) >= min_valid) %>%
  arrange(desc(Latent_Magnitude))

increased_set <- df_sd %>%
  filter(Direction == "INCREASED") %>%
  filter(rowSums(!is.na(select(., all_of(fc_increased_cols)))) >= min_valid) %>%
  arrange(desc(Latent_Magnitude))

decreased_top <- head(decreased_set, n_take)
increased_top <- head(increased_set, n_take)

# 8) Export selected tables

write.xlsx(
  decreased_top,
  file.path(out_tables, paste0("Selected_DECREASED_top", nrow(decreased_top), "_", today, ".xlsx"))
)

write.xlsx(
  increased_top,
  file.path(out_tables, paste0("Selected_INCREASED_top", nrow(increased_top), "_", today, ".xlsx"))
)

# 9) Split into groups of 10 (for heatmaps per direction)
chunk10 <- function(d) {
  if (nrow(d) == 0) return(list())
  idx <- rep(1:ceiling(nrow(d) / 10), each = 10)[seq_len(nrow(d))]
  split(d, idx)
}
decreased_chunks <- chunk10(decreased_top)
increased_chunks <- chunk10(increased_top)

# 10) Heatmap function (robust, consistent, styled)

make_heatmap <- function(dsub, cols_fc, title_prefix, file_stub, out_dir_fig) {

  if (nrow(dsub) == 0) return(invisible(NULL))

  # Expect coverage of 4 comparisons (female & male)
  expected_groups <- c("FC_FCvsFW", "FC_FWRvsFW", "FC_MCvsMW", "FC_MWRvsMW")

  # Keep only FC columns that match expected groups AND belong to this direction (cols_fc)
  expected_cols <- grep(paste(expected_groups, collapse = "|"), cols_fc, value = TRUE)

  # Ensure missing expected columns exist (as NA)
  for (ec in expected_cols) {
    if (!ec %in% names(dsub)) dsub[[ec]] <- NA
  }

  cols_fc_use <- intersect(expected_cols, names(dsub))
  if (length(cols_fc_use) == 0) return(invisible(NULL))

  mat_fc <- as.matrix(dsub[, cols_fc_use, drop = FALSE])
  rownames(mat_fc) <- dsub$Accession

  # Treat 0 as missing (consistent with your original logic)
  mat_fc[mat_fc == 0] <- NA

  # log2 transform
  mat_log2 <- log2(mat_fc)

  # Drop empty rows/cols
  keep_rows <- rowSums(is.na(mat_log2)) < ncol(mat_log2)
  keep_cols <- colSums(is.na(mat_log2)) < nrow(mat_log2)
  mat_log2 <- mat_log2[keep_rows, keep_cols, drop = FALSE]
  dsub_keep <- dsub[keep_rows, , drop = FALSE]

  if (nrow(mat_log2) == 0 || ncol(mat_log2) == 0) return(invisible(NULL))

  # Cap extremes for display
  cap <- 2
  mat_disp <- mat_log2
  mat_disp[mat_disp >  cap] <-  cap
  mat_disp[mat_disp < -cap] <- -cap

  # Median imputation for clustering only
  col_meds <- apply(mat_log2, 2, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) 0 else median(x, na.rm = TRUE)
  })

  mat_clust <- mat_log2
  for (j in seq_len(ncol(mat_clust))) {
    idx <- which(!is.finite(mat_clust[, j]) | is.na(mat_clust[, j]))
    if (length(idx)) mat_clust[idx, j] <- col_meds[j]
  }

  cl_rows <- if (nrow(mat_clust) >= 2) hclust(dist(mat_clust)) else FALSE
  cl_cols <- if (ncol(mat_clust) >= 2) hclust(dist(t(mat_clust))) else FALSE

  # Colors & annotations
  col_fun <- circlize::colorRamp2(c(-cap, 0, cap), c("#2b6cb0", "white", "#e34a33"))

  ha_row <- rowAnnotation(
    LatentMag = anno_barplot(
      scale(dsub_keep$Latent_Magnitude),
      gp = gpar(fill = "#4a5568", col = NA),
      border = FALSE, width = unit(2.5, "cm"),
      axis_param = list(gp = gpar(cex = 0.6))
    ),
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontface = "bold")
  )

  col_group <- ifelse(grepl("^FC_FC|^FC_FWR", colnames(mat_disp)), "Female",
                      ifelse(grepl("^FC_MC|^FC_MWR", colnames(mat_disp)), "Male", "Other"))

  ca <- HeatmapAnnotation(
    Group = col_group,
    col = list(Group = c(Female = "#2563eb", Male = "#ef4444", Other = "#737373")),
    annotation_name_gp = gpar(fontface = "bold")
  )

  ht <- Heatmap(
    mat_disp,
    name = "log2(FC)",
    col = col_fun,
    na_col = "#e5e7eb",
    top_annotation = ca,
    right_annotation = ha_row,
    cluster_rows = cl_rows,
    cluster_columns = cl_cols,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(cex = 0.65),
    column_names_gp = gpar(cex = 0.75, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-2, -1, 0, 1, 2),
      title_gp = gpar(fontface = "bold"),
      labels_gp = gpar(cex = 0.8)
    ),
    column_title = paste0(title_prefix, " (n=", nrow(mat_disp), ")"),
    column_title_gp = gpar(fontface = "bold", cex = 1.0)
  )

  out_file <- file.path(out_dir_fig, paste0(file_stub, "_", today, ".png"))
  png(filename = out_file, width = 1800, height = 1300, res = 220)
  draw(ht, padding = unit(c(6, 6, 6, 6), "mm"))
  dev.off()

  message("Saved: ", out_file)
}

# 11) Generate 6 heatmaps (3 DECREASED + 3 INCREASED)

labels <- c("1–10", "11–20", "21–30")

for (i in seq_along(decreased_chunks)) {
  make_heatmap(
    decreased_chunks[[i]],
    fc_decreased_cols,
    title_prefix = paste0("DECREASED protein abundance — Group ", i, " (Top ", labels[i], ")"),
    file_stub = paste0("Heatmap_DECREASED_Group", i),
    out_dir_fig = out_figs
  )
}

for (i in seq_along(increased_chunks)) {
  make_heatmap(
    increased_chunks[[i]],
    fc_increased_cols,
    title_prefix = paste0("INCREASED protein abundance — Group ", i, " (Top ", labels[i], ")"),
    file_stub = paste0("Heatmap_INCREASED_Group", i),
    out_dir_fig = out_figs
  )
}

# 12) Summary
cat("\n Latent screening + heatmap generation complete!\n")
cat("DECREASED proteins selected:", nrow(decreased_top), "→ heatmaps:", length(decreased_chunks), "\n")
cat("INCREASED proteins selected:", nrow(increased_top), "→ heatmaps:", length(increased_chunks), "\n")
cat("Tables saved to:", out_tables, "\n")
cat("Figures saved to:", out_figs, "\n")
