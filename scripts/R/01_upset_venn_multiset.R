# ==========================================================
# UpSet + Venn Plot Generator for protein abundance Sets
# - UpSet (ComplexUpset) with vertical labels
# - Split Venns + unified 7-set Venn (ggVennDiagram)
# Manuscript: Deep learning-assisted proteomic dissection reveals sex biased shared proteomic patterns in Populus deltoides under waterlogging stress and recovery  
# Author: Dr. El Hadji Malick Cisse
# ==========================================================
  library(readxl)
  library(dplyr)
  library(ComplexUpset)
  library(ggplot2)
  library(ggVennDiagram)
  library(patchwork)

# Paths
up_file   <- file.path("data", "Increased.xlsx")
down_file <- file.path("data", "Decreased.xlsx")

# Output directory
out_dir <- file.path("results", "figures", "upset_venn")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: Read and clean Excel gene sets 
read_and_clean <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste0(
      "Input file not found: ", file_path, "\n",
      "Place the file in the expected location or edit the path."
    ))
  }

  df <- readxl::read_excel(file_path)

  # Remove duplicated column names, then make names unique
  df <- df[, !duplicated(names(df)), drop = FALSE]
  names(df) <- make.unique(names(df))

  # Each column becomes a unique vector of proteins IDs
  lapply(df, function(x) unique(na.omit(as.character(x))))
}

# Binary membership matrix
make_binary_df <- function(list_data) {
  all_proteins <- unique(unlist(list_data))
  df_bin <- as.data.frame(sapply(list_data, function(set) all_proteins %in% set))
  rownames(df_bin) <- all_proteins
  df_bin
}

# Biological order
bio_order <- c("FW vs FC", "MW vs MC", "MW vs FW", "MWR vs MW",
               "MWR vs FWR", "FWR vs FW", "MC vs FC")

# Read data
up_sets   <- read_and_clean(up_file)
down_sets <- read_and_clean(down_file)

# Colors
cols_up <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
             "#e6ab02","#a6761d","#666666")
cols_down <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3",
               "#fdb462","#b3de69","#fccde5")

# Build binary matrices and order columns
Increased_bin <- make_binary_df(up_sets)
Decreased_bin <- make_binary_df(down_sets)

Increased_bin <- up_bin[, intersect(bio_order, colnames(up_bin)), drop = FALSE]
Decreased_bin <- down_bin[, intersect(bio_order, colnames(down_bin)), drop = FALSE]

# UpSet plot function
plot_upset <- function(df_bin, title_text, fill_cols, out_png) {
  png(out_png, width = 3000, height = 2000, res = 300)

  p <- ComplexUpset::upset(
    df_bin,
    names(df_bin),
    height_ratio = c(0.65, 0.35),
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        text = list(size = 2.8, color = "black", face = "bold", vjust = -0.5)
      )
    )
  ) +
    labs(title = title_text, y = "Intersection size", x = "") +
    scale_fill_manual(values = fill_cols) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(
        size = 6.5, color = "gray25",
        angle = 90, hjust = 1, vjust = 0.5,
        margin = margin(t = 6)
      ),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )

  print(p)
  dev.off()
}

# Save UpSet plots                                
plot_upset(Increased_bin, "Increased proteins", cols_up,
           file.path(out_dir, "Increased_UpSet.png"))

plot_upset(Decreased_bin, "Decreased proteins", cols_down,
           file.path(out_dir, "Decreased_UpSet.png"))

# Helpers for Venn sets
make_unique_sets <- function(list_data) {
  list_data[!grepl("\\.\\.\\.\\d+$", names(list_data))]
}

split_sets <- function(list_data) {
  list_data <- make_unique_sets(list_data)
  n <- length(list_data)
  mid <- ceiling(n / 2)
  list(Group1 = list_data[1:mid], Group2 = list_data[(mid + 1):n])
}

increased_groups <- split_sets(increased_sets)
decreased_groups <- split_sets(decreased_sets)

# Split Venns (Increased)
png(file.path(out_dir, "Increased_Venn_split.png"), width = 2500, height = 1200, res = 300)
p1 <- ggVennDiagram(increased_groups$Group1, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#e7298a") +
  ggtitle("Increased (Group 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
p2 <- ggVennDiagram(Increased_groups$Group2, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#1b9e77") +
  ggtitle("Increased (Group 2)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
print(p1 + p2)
dev.off()

# Split Venns (Decreased)
png(file.path(out_dir, "Decreased_Venn_split.png"), width = 2500, height = 1200, res = 300)
p3 <- ggVennDiagram(Decreased_groups$Group1, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#fb8072") +
  ggtitle("Decreased (Group 1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
p4 <- ggVennDiagram(down_groups$Group2, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#80b1d3") +
  ggtitle("Decreased (Group 2)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
print(p3 + p4)
dev.off()

# Unified 7-set Venns
Increased7 <- make_unique_sets(Increased_sets)
Decreased7 <- make_unique_sets(Decreased_sets)

png(file.path(out_dir, "Increased_Venn7.png"), width = 2500, height = 2000, res = 300)
ggVennDiagram(up7, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#66a61e") +
  ggtitle("Increased (All comparisons)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
dev.off()

png(file.path(out_dir, "Decreased_Venn7.png"), width = 2500, height = 2000, res = 300)
ggVennDiagram(down7, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#80b1d3") +
  ggtitle("Decreased (All comparisons)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(color = "black", face = "bold"))
dev.off()

message(" UpSet + Venn plots saved to: ", out_dir)
