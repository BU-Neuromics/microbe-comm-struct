#!/usr/bin/env Rscript
# =============================================================================
# Script: 01_aitchison_pca.R
#
# Description:
#   Aitchison PCA (CLR-transformed PCA) for microbial community structure
#   analysis. Zeros are replaced using zCompositions::cmultRepl (CZM method)
#   before CLR transformation. PCA is run on the samples x taxa CLR matrix.
#
# Inputs:
#   1. counts_csv   - Path to counts matrix CSV (rows = taxa, cols = samples)
#   2. metadata_csv - Path to metadata CSV (rows = samples, must have 'group')
#
# Outputs (saved to same directory as counts_csv):
#   - aitchison_scores.pdf  - PCA scores plot colored by group
#   - aitchison_loadings.pdf - Top 30 taxa loadings plot
#
# Stdout:
#   - Variance explained per PC (first 10 PCs)
#   - Pearson correlations of PC1/PC2 scores with library size
#
# Dependencies:
#   zCompositions, ggplot2, ggrepel
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(zCompositions)
  library(ggplot2)
  library(ggrepel)
})

# --- Parse arguments ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 01_aitchison_pca.R <counts_csv> <metadata_csv>\n",
      file = stderr())
  quit(status = 1)
}

counts_path   <- args[1]
metadata_path <- args[2]
out_dir       <- dirname(counts_path)

# --- Load data ----------------------------------------------------------------
message("Loading counts matrix from: ", counts_path)
counts_raw <- read.csv(counts_path, row.names = 1, check.names = FALSE)
counts_raw <- as.matrix(counts_raw)
mode(counts_raw) <- "numeric"
# Input is samples x taxa; transpose to taxa x samples for internal use
counts_raw <- t(counts_raw)

message("Loading metadata from: ", metadata_path)
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

# Align samples: columns of counts (now taxa x samples) match rows of metadata
common_samples <- intersect(colnames(counts_raw), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No samples in common between counts rows and metadata rows.",
       call. = FALSE)
}
message("Samples found in both counts and metadata: ", length(common_samples))
counts_raw <- counts_raw[, common_samples, drop = FALSE]
metadata   <- metadata[common_samples, , drop = FALSE]

if (!"group" %in% colnames(metadata)) {
  stop("Metadata must have a column named 'group'.", call. = FALSE)
}

# --- Library sizes ------------------------------------------------------------
lib_sizes <- colSums(counts_raw)
message("Library sizes range: [", round(min(lib_sizes), 2), ", ",
        round(max(lib_sizes), 2), "]")

# --- Zero replacement with cmultRepl (CZM) ------------------------------------
# cmultRepl expects samples as rows, taxa as columns
message("Replacing zeros using zCompositions::cmultRepl (CZM)...")
counts_t <- t(counts_raw)  # samples x taxa

# cmultRepl requires a numeric matrix; fractional counts are supported
counts_nozero <- cmultRepl(counts_t, method = "CZM", output = "p-counts",
                            suppress.print = TRUE)
message("Zero replacement complete. Matrix dimensions: ",
        nrow(counts_nozero), " samples x ", ncol(counts_nozero), " taxa")

# --- CLR transformation -------------------------------------------------------
message("Applying CLR transformation...")
# CLR: log(x) - mean(log(x)) for each sample (row)
log_counts <- log(counts_nozero)
row_log_means <- rowMeans(log_counts)
clr_matrix <- log_counts - row_log_means  # broadcast: subtract per-row mean

# --- PCA ----------------------------------------------------------------------
message("Running PCA on CLR-transformed matrix (", nrow(clr_matrix),
        " samples x ", ncol(clr_matrix), " taxa)...")
pca_result <- prcomp(clr_matrix, center = TRUE, scale. = FALSE)

# Variance explained
var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2)
n_pcs_to_show <- min(10, length(var_explained))

cat("\n=== Variance Explained (first 10 PCs) ===\n")
ve_table <- data.frame(
  PC = paste0("PC", seq_len(n_pcs_to_show)),
  Variance_Explained = round(var_explained[seq_len(n_pcs_to_show)], 4),
  Cumulative = round(cumsum(var_explained)[seq_len(n_pcs_to_show)], 4)
)
print(ve_table, row.names = FALSE)

# Pearson correlations of PC1, PC2 with library size
pc1_scores <- pca_result$x[, 1]
pc2_scores <- pca_result$x[, 2]
cor_pc1_libsize <- cor(pc1_scores, lib_sizes, method = "pearson")
cor_pc2_libsize <- cor(pc2_scores, lib_sizes, method = "pearson")

cat("\n=== Pearson Correlation of PC Scores with Library Size ===\n")
cat("PC1 vs library size: r =", round(cor_pc1_libsize, 4), "\n")
cat("PC2 vs library size: r =", round(cor_pc2_libsize, 4), "\n\n")

# --- Build scores data frame --------------------------------------------------
scores_df <- data.frame(
  sample_id  = rownames(pca_result$x),
  PC1        = pca_result$x[, 1],
  PC2        = pca_result$x[, 2],
  group      = metadata[rownames(pca_result$x), "group"],
  lib_size   = lib_sizes[rownames(pca_result$x)],
  stringsAsFactors = FALSE
)

pc1_pct <- round(var_explained[1] * 100, 1)
pc2_pct <- round(var_explained[2] * 100, 1)

# --- Scores plot --------------------------------------------------------------
message("Generating scores plot...")
scores_plot <- ggplot(scores_df, aes(x = PC1, y = PC2,
                                     color = group,
                                     size  = lib_size,
                                     label = sample_id)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_size_continuous(name = "Library Size",
                        labels = scales::comma) +
  labs(
    title  = "Aitchison PCA — Scores",
    x      = paste0("PC1 (", pc1_pct, "% variance explained)"),
    y      = paste0("PC2 (", pc2_pct, "% variance explained)"),
    color  = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

scores_out <- file.path(out_dir, "aitchison_scores.pdf")
ggsave(scores_out, scores_plot, width = 9, height = 7)
message("Scores plot saved to: ", scores_out)

# --- Loadings plot ------------------------------------------------------------
message("Generating loadings plot...")
loadings_mat <- pca_result$rotation  # taxa x PCs

# Top 30 taxa by sum of squared loadings on PC1 + PC2
ss_loadings <- loadings_mat[, 1]^2 + loadings_mat[, 2]^2
top30_idx    <- order(ss_loadings, decreasing = TRUE)[seq_len(min(30, nrow(loadings_mat)))]
top30_taxa   <- rownames(loadings_mat)[top30_idx]

loadings_df <- data.frame(
  taxon  = top30_taxa,
  PC1    = loadings_mat[top30_taxa, 1],
  PC2    = loadings_mat[top30_taxa, 2],
  ss     = ss_loadings[top30_taxa],
  stringsAsFactors = FALSE
)
# Order by ss for plotting
loadings_df <- loadings_df[order(loadings_df$ss, decreasing = TRUE), ]
loadings_df$taxon <- factor(loadings_df$taxon, levels = rev(loadings_df$taxon))

loadings_plot <- ggplot(loadings_df, aes(x = taxon, y = ss)) +
  geom_point(aes(color = PC1), size = 3) +
  geom_segment(aes(xend = taxon, yend = 0), color = "grey50", linewidth = 0.4) +
  coord_flip() +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato",
                        midpoint = 0, name = "PC1 Loading") +
  labs(
    title = "Aitchison PCA — Top 30 Taxa Loadings",
    x     = "Taxon",
    y     = "Sum of Squared Loadings (PC1 + PC2)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7)
  )

loadings_out <- file.path(out_dir, "aitchison_loadings.pdf")
ggsave(loadings_out, loadings_plot, width = 9, height = 10)
message("Loadings plot saved to: ", loadings_out)

message("Script 01_aitchison_pca.R complete.")
