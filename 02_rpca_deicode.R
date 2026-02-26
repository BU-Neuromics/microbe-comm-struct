#!/usr/bin/env Rscript
# =============================================================================
# Script: 02_rpca_deicode.R
#
# Description:
#   Robust PCA (RPCA / DEICODE approach) for microbial community analysis.
#   Uses robust CLR (rCLR) — CLR computed using only non-zero entries per
#   sample — followed by matrix completion via softImpute to fill zero
#   positions, then PCA via rsvd for efficiency on large taxa sets.
#   Rank is selected by cross-validated held-out reconstruction error.
#
# Inputs:
#   1. counts_csv   - Path to counts matrix CSV (rows = taxa, cols = samples)
#   2. metadata_csv - Path to metadata CSV (rows = samples, must have 'group')
#
# Outputs (saved to same directory as counts_csv):
#   - rpca_scores.pdf   - PCA scores plot colored by group
#   - rpca_loadings.pdf - Top 30 taxa loadings plot
#
# Stdout:
#   - Chosen rank and held-out reconstruction error
#
# Dependencies:
#   softImpute, rsvd, ggplot2, ggrepel
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(softImpute)
  library(rsvd)
  library(ggplot2)
  library(ggrepel)
})

# --- Parse arguments ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 02_rpca_deicode.R <counts_csv> <metadata_csv>\n",
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

# --- rCLR transformation ------------------------------------------------------
# rCLR: for each sample (column), compute CLR using only non-zero entries.
# Zero positions are set to NA (to be filled by matrix completion).
message("Applying robust CLR (rCLR) transformation...")

# Work in samples x taxa orientation (transpose)
counts_t <- t(counts_raw)  # samples x taxa
n_samples <- nrow(counts_t)
n_taxa    <- ncol(counts_t)

rclr_matrix <- matrix(NA_real_, nrow = n_samples, ncol = n_taxa)
rownames(rclr_matrix) <- rownames(counts_t)
colnames(rclr_matrix) <- colnames(counts_t)

for (i in seq_len(n_samples)) {
  row_i <- counts_t[i, ]
  nonzero_idx <- which(row_i > 0)
  if (length(nonzero_idx) == 0) {
    # All zeros: leave as NA
    next
  }
  log_nz <- log(row_i[nonzero_idx])
  gm_nz  <- mean(log_nz)  # log geometric mean of non-zero entries
  rclr_matrix[i, nonzero_idx] <- log_nz - gm_nz
  # Zero positions remain NA
}

# Record positions of zero entries (will be used for cross-validation)
zero_mask <- is.na(rclr_matrix)
n_zeros   <- sum(zero_mask)
n_nonzero <- sum(!zero_mask)
message("rCLR complete. Non-zero entries: ", n_nonzero,
        ", Zero (NA) positions: ", n_zeros)

# --- Rank selection by cross-validation on non-zero entries -------------------
message("Selecting rank via cross-validated held-out reconstruction error...")
message("  (20% random holdout of non-zero entries, testing ranks 1 to ",
        floor(min(n_samples, n_taxa) / 2), ")")

max_rank <- floor(min(n_samples, n_taxa) / 2)
max_rank <- min(max_rank, 20)  # cap at 20 for computational feasibility

# Identify non-zero positions as linear indices
nonzero_positions <- which(!zero_mask)
holdout_n         <- floor(0.2 * length(nonzero_positions))
set.seed(42)
holdout_idx       <- sample(nonzero_positions, holdout_n)

# Build training matrix: mask holdout positions as NA
rclr_train <- rclr_matrix
rclr_train[holdout_idx] <- NA

holdout_true <- rclr_matrix[holdout_idx]

cv_errors <- numeric(max_rank)

for (k in seq_len(max_rank)) {
  fit_k <- tryCatch(
    softImpute(rclr_train, rank.max = k, lambda = 0, type = "svd",
               trace.it = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit_k)) {
    cv_errors[k] <- Inf
    next
  }
  completed_k  <- complete(rclr_train, fit_k)
  holdout_pred <- completed_k[holdout_idx]
  cv_errors[k] <- sqrt(mean((holdout_true - holdout_pred)^2))
}

best_rank <- which.min(cv_errors)
best_rmse <- cv_errors[best_rank]

cat("\n=== Rank Selection Results ===\n")
cat("Chosen rank:", best_rank, "\n")
cat("Held-out reconstruction RMSE:", round(best_rmse, 6), "\n\n")

# Print full CV error table
cv_df <- data.frame(rank = seq_len(max_rank), cv_rmse = round(cv_errors, 6))
print(cv_df, row.names = FALSE)
cat("\n")

# --- Matrix completion at chosen rank -----------------------------------------
message("Running matrix completion at rank ", best_rank, "...")
fit_final    <- softImpute(rclr_matrix, rank.max = best_rank, lambda = 0,
                           type = "svd", trace.it = FALSE)
rclr_completed <- complete(rclr_matrix, fit_final)
message("Matrix completion done. Dimensions: ", nrow(rclr_completed),
        " x ", ncol(rclr_completed))

# --- PCA via rsvd -------------------------------------------------------------
message("Running randomized SVD (rsvd) on completed matrix...")
k_svd      <- min(best_rank + 2, min(n_samples, n_taxa) - 1)
svd_result <- rsvd(rclr_completed, k = k_svd)

# Compute scores and loadings
# rsvd returns: u (samples x k), d (k), v (taxa x k)
scores_mat   <- svd_result$u %*% diag(svd_result$d)
loadings_mat <- svd_result$v
rownames(scores_mat)   <- rownames(rclr_completed)
rownames(loadings_mat) <- colnames(rclr_completed)
colnames(scores_mat)   <- paste0("PC", seq_len(ncol(scores_mat)))
colnames(loadings_mat) <- paste0("PC", seq_len(ncol(loadings_mat)))

# Approximate variance explained from singular values
total_ss     <- sum(rclr_completed^2)
var_explained <- svd_result$d^2 / total_ss

n_pcs_show <- min(ncol(scores_mat), 10)
cat("=== Approximate Variance Explained (RPCA) ===\n")
ve_table <- data.frame(
  PC = paste0("PC", seq_len(n_pcs_show)),
  Variance_Explained = round(var_explained[seq_len(n_pcs_show)], 4),
  Cumulative = round(cumsum(var_explained)[seq_len(n_pcs_show)], 4)
)
print(ve_table, row.names = FALSE)
cat("\n")

# --- Scores plot --------------------------------------------------------------
message("Generating scores plot...")
pc1_pct <- round(var_explained[1] * 100, 1)
pc2_pct <- round(var_explained[2] * 100, 1)

scores_df <- data.frame(
  sample_id = rownames(scores_mat),
  PC1       = scores_mat[, 1],
  PC2       = scores_mat[, 2],
  group     = metadata[rownames(scores_mat), "group"],
  lib_size  = lib_sizes[rownames(scores_mat)],
  stringsAsFactors = FALSE
)

scores_plot <- ggplot(scores_df, aes(x = PC1, y = PC2,
                                     color = group,
                                     size  = lib_size,
                                     label = sample_id)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_size_continuous(name = "Library Size",
                        labels = scales::comma) +
  labs(
    title = "Robust PCA (DEICODE/rCLR) — Scores",
    x     = paste0("PC1 (", pc1_pct, "% var. exp.)"),
    y     = paste0("PC2 (", pc2_pct, "% var. exp.)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "right"
  )

scores_out <- file.path(out_dir, "rpca_scores.pdf")
ggsave(scores_out, scores_plot, width = 9, height = 7)
message("Scores plot saved to: ", scores_out)

# --- Loadings plot ------------------------------------------------------------
message("Generating loadings plot...")
ss_loadings <- loadings_mat[, 1]^2 + loadings_mat[, 2]^2
top30_idx    <- order(ss_loadings, decreasing = TRUE)[seq_len(min(30, nrow(loadings_mat)))]
top30_taxa   <- rownames(loadings_mat)[top30_idx]

loadings_df <- data.frame(
  taxon = top30_taxa,
  PC1   = loadings_mat[top30_taxa, 1],
  PC2   = loadings_mat[top30_taxa, 2],
  ss    = ss_loadings[top30_taxa],
  stringsAsFactors = FALSE
)
loadings_df <- loadings_df[order(loadings_df$ss, decreasing = TRUE), ]
loadings_df$taxon <- factor(loadings_df$taxon, levels = rev(loadings_df$taxon))

loadings_plot <- ggplot(loadings_df, aes(x = taxon, y = ss)) +
  geom_point(aes(color = PC1), size = 3) +
  geom_segment(aes(xend = taxon, yend = 0), color = "grey50", linewidth = 0.4) +
  coord_flip() +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato",
                        midpoint = 0, name = "PC1 Loading") +
  labs(
    title = "Robust PCA — Top 30 Taxa Loadings",
    x     = "Taxon",
    y     = "Sum of Squared Loadings (PC1 + PC2)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 7)
  )

loadings_out <- file.path(out_dir, "rpca_loadings.pdf")
ggsave(loadings_out, loadings_plot, width = 9, height = 10)
message("Loadings plot saved to: ", loadings_out)

message("Script 02_rpca_deicode.R complete.")
