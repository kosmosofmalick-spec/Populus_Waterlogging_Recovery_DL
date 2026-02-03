# ==========================================================
# Sex dimorphism (scale-free) — ggplot figures
# Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
# Reads CSVs produced by:
# scripts/Python/04_sex_dimorphism_scale_free_contrasts.py
# Outputs:
# - results/figures/sex_dimorphism/*.png
# Author: Dr. El Hadji Malick Cisse
# ==========================================================

  library(tidyverse)
  library(ggrepel)
  library(scales)
  library(readr)
  library(patchwork)

base_tables <- file.path("results", "tables", "sex_dimorphism")
out_figs <- file.path("results", "figures", "sex_dimorphism")
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

p_ <- function(x) file.path(base_tables, x)
f_ <- function(x) file.path(out_figs, x)

# 1) Load data
contrib_all <- read_csv(p_("Z_contributions_all_contrasts.csv"), show_col_types = FALSE)
sex_diffs   <- read_csv(p_("Z_sex_diffs_by_treatment.csv"), show_col_types = FALSE)
means_clean <- read_csv(p_("traits_means_clean.csv"), show_col_types = FALSE)

# 2) Figure 1 — Faceted contribution bars
contrib_ordered <- contrib_all %>%
  group_by(contrast) %>%
  mutate(Trait = fct_reorder(Trait, percent_contribution, .fun = identity)) %>%
  ungroup()

fig1 <- contrib_ordered %>%
  ggplot(aes(x = percent_contribution, y = Trait)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percent_contribution)),
            hjust = -0.05, size = 3.2, family = "sans") +
  facet_wrap(~contrast, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Traits characterizing sex dimorphism and recovery (standardized, scale-free)",
    x = "Percent contribution (mean |Δ|, z-scored)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(f_("Z_contributions_faceted_barplots.png"), fig1, width = 11, height = 7, dpi = 300)

# 3) Figure 2 — PCA on standardized means
id_cols <- c("Treatments", "Sex")
trait_cols <- setdiff(names(means_clean), id_cols)

means_z <- means_clean %>%
  mutate(across(all_of(trait_cols), ~ (. - mean(.)) / sd(.)))

pca <- prcomp(means_z %>% select(all_of(trait_cols)), center = FALSE, scale. = FALSE)

scores <- as_tibble(pca$x) %>%
  bind_cols(means_z %>% select(all_of(id_cols)))

ve <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", percent(ve[1], accuracy = 0.1), ")")
pc2_lab <- paste0("PC2 (", percent(ve[2], accuracy = 0.1), ")")

fig2 <- ggplot(scores, aes(PC1, PC2, shape = Sex, label = Treatments)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, min.segment.length = 0, family = "sans") +
  labs(
    title = "PCA on standardized trait means (Treatment × Sex)",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(f_("Z_PCA_sex_treatments.png"), fig2, width = 7, height = 5, dpi = 300)

# 4) Figure 3 — Heatmap of standardized Male − Female diffs per treatment
sex_diffs_long <- sex_diffs %>%
  pivot_longer(-Treatment, names_to = "Trait", values_to = "z_diff")

fig3 <- ggplot(sex_diffs_long, aes(Treatment, Trait, fill = z_diff)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, guide = guide_colorbar(barheight = unit(60, "pt"))) +
  labs(
    title = "Standardized sex differences per treatment (Male − Female)",
    x = NULL, y = NULL, fill = "z-diff"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())

ggsave(f_("Z_heatmap_sex_diffs.png"), fig3, width = 7, height = 5, dpi = 300)

# 5) Figure 4 — Sex Dimorphism Index (SDI) per treatment
sdi_tbl <- sex_diffs %>%
  mutate(SDI_abs_mean = rowMeans(across(-Treatment, ~ abs(.)))) %>%
  select(Treatment, SDI_abs_mean)

fig4 <- sdi_tbl %>%
  ggplot(aes(x = SDI_abs_mean, y = fct_reorder(Treatment, SDI_abs_mean))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(SDI_abs_mean, 2)),
            hjust = -0.05, size = 3.3, family = "sans") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Sex Dimorphism Index (SDI, mean |z-diff|) by treatment",
    x = "Mean |Male − Female| (z-score units)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(f_("Z_SDI_by_treatment.png"), fig4, width = 7, height = 4.5, dpi = 300)

# 6) Combo panel
combo <- fig1 / fig4 + plot_layout(heights = c(2, 1))
ggsave(f_("Z_contributions_plus_SDI_combo.png"), combo, width = 11, height = 10, dpi = 300)

cat("Saved figures in: ", out_figs, "\n", sep = "")
