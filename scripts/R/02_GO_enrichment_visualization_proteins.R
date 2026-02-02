# ==========================================================
# GO Enrichment Visualization (Bar, Bubble, Treemap)
# Proteomics / protein abundance
# Manuscript: Deep learning-assisted proteomic dissection reveals sex biased and shared proteomic patterns in Populus deltoides under waterlogging stress and recovery
# Author: Dr. El Hadji Malick Cisse
# ==========================================================

  library(tidyverse)
  library(readxl)
  library(ggrepel)
  library(treemapify)

# 1) Input/Output

input_path <- file.path("data", "GO_enrichment.xlsx")

out_dir <- file.path("results", "figures", "go_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop(paste0(
    "Input file not found: ", input_path, "\n",
    "This GitHub repo is code-only; place GO enrichment Excel file in /data\n",
    "or edit `input_path` to your local location."
  ))
}

# 2) Read Excel

df <- readxl::read_excel(input_path)

# 3) Clean and preprocess

# NOTE: column names below must match your Excel export.
# If the file uses slightly different headers, edit the rename() mapping.
df_clean <- df %>%
  rename(
    GO_ID = `Gene Ontology ID`,
    Term = `Gene Ontology term`,
    ClusterFreq = `Cluster frequency`,
    Background = `Protein frequency of use`,
    Pvalue = `P-value`,
    Proteins = `Genes annotated to the term`
  ) %>%
  mutate(
    # Extract numeric count from "38 of 184 in the list"
    ClusterCount = as.numeric(sub(" of.*", "", as.character(ClusterFreq))),
    ProteinCount = ifelse(is.na(ClusterCount), 0, ClusterCount),
    logP = -log10(as.numeric(Pvalue))
  ) %>%
  # avoid duplicated terms (keep most significant instance)
  arrange(Pvalue) %>%
  distinct(Term, .keep_all = TRUE)

# 
# 4) Select top enriched terms
# 
df_top <- df_clean %>%
  filter(!is.na(Pvalue), is.finite(logP)) %>%
  arrange(Pvalue) %>%
  slice_head(n = 15) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term))))

# ==========================================================
# BAR PLOT (Top Enriched GO Terms)
# ==========================================================
bar_plot <- ggplot(df_top, aes(x = logP, y = Term, fill = logP)) +
  geom_col(width = 0.7, color = "white") +
  labs(
    x = expression(-log[10](P[value])),
    y = "GO term",
    title = "Top enriched GO terms (proteins)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "GO_Barplot.pdf"), bar_plot, width = 7, height = 5)
ggsave(file.path(out_dir, "GO_Barplot.png"), bar_plot, width = 7, height = 5, dpi = 600)

# ==========================================================
# BUBBLE PLOT (Significance vs Protein Count)
# ==========================================================
bubble_plot <- ggplot(df_top, aes(x = logP, y = Term, size = ProteinCount, color = logP)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(3, 10)) +
  labs(
    x = expression(-log[10](P[value])),
    y = "GO term",
    color = expression(-log[10](P[value])),
    size = "Protein count",
    title = "GO enrichment bubble plot (proteins)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key.size = unit(0.6, "cm")
  )

ggsave(file.path(out_dir, "GO_BubblePlot.pdf"), bubble_plot, width = 7, height = 5)
ggsave(file.path(out_dir, "GO_BubblePlot.png"), bubble_plot, width = 7, height = 5, dpi = 600)

# ==========================================================
# TREEMAP (Hierarchy of Enriched Terms)
# ==========================================================
treemap_data <- df_top %>%
  mutate(Label = paste0(Term, "\n(", ProteinCount, " proteins)"))

treemap_plot <- ggplot(treemap_data,
                       aes(area = ProteinCount, fill = logP, label = Label)) +
  geom_treemap(color = "white") +
  geom_treemap_text(
    fontface = "italic", colour = "white",
    place = "centre", reflow = TRUE, min.size = 6
  ) +
  labs(
    title = "GO enrichment treemap (proteins)",
    fill = expression(-log[10](P[value]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(file.path(out_dir, "GO_Treemap.pdf"), treemap_plot, width = 7, height = 6)
ggsave(file.path(out_dir, "GO_Treemap.png"), treemap_plot, width = 7, height = 6, dpi = 600)

message(" GO enrichment figures saved to: ", out_dir)
