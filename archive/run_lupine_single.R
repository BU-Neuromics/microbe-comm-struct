#!/usr/bin/env Rscript

# Script to run LUPINE (single time point mode) on a CSV count file.
#
# LUPINE uses partial least squares regression to infer microbial co-occurrence
# networks.  The single time point mode treats the data as one snapshot and
# computes pairwise partial correlations after removing shared low-dimensional
# structure (PCA/ICA/RPCA components).
#
# Outputs mirror run_netcomi_method.R so that compare_association_methods.R can
# consume them without modification.

# ── Package checks ────────────────────────────────────────────────────────────

if (!requireNamespace("optparse", quietly = TRUE)) {
  message("Installing optparse...")
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}

if (!requireNamespace("LUPINE", quietly = TRUE)) {
  message("Installing LUPINE from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org/")
  }
  devtools::install_github("https://github.com/SarithaKodikara/LUPINE")
}

library(LUPINE)
library(optparse)

# ── CLI options ───────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to CSV file: samples as rows, taxa as columns [required]",
              metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "lupine_output",
              help = "Prefix for output files [default: %default]",
              metavar = "PREFIX"),
  make_option(c("-p", "--prevalence"), type = "double", default = 0.5,
              help = "Minimum prevalence threshold (proportion of samples) [default: %default]",
              metavar = "PROP"),
  make_option("--minsamp", type = "integer", default = 1000,
              help = "Minimum total read count per sample; samples below are dropped [default: %default]",
              metavar = "N"),
  make_option("--maxtaxa", type = "integer", default = 10000,
              help = "Keep at most this many taxa (top by mean abundance after prevalence filter) [default: %default]",
              metavar = "N"),
  make_option("--dimred", type = "character", default = "pca",
              help = "Dimensionality reduction method: pca, ica, rpca [default: %default]",
              metavar = "METHOD"),
  make_option("--ncomp", type = "integer", default = 1,
              help = "Number of components for dimensionality reduction [default: %default]",
              metavar = "N"),
  make_option("--cutoff", type = "double", default = 0.05,
              help = "p-value cutoff for declaring a significant edge [default: %default]",
              metavar = "PVAL"),
  make_option("--transformed", type = "logical", default = FALSE,
              help = "Set TRUE if data are already CLR- or log-transformed (skips library-size offset) [default: %default]",
              metavar = "BOOL"),
  make_option("--seed", type = "integer", default = 123456,
              help = "Random seed for reproducibility [default: %default]",
              metavar = "SEED")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nRun LUPINE single time point network inference on taxonomic count data.\n\nThe input CSV must have samples as rows and taxa as columns (same format\nas used by run_netcomi_method.R and run_spring_from_csv.R).",
  epilogue = "Example:\n  Rscript run_lupine_single.R -i counts.csv -o lupine_results --cutoff 0.05"
)
opt <- parse_args(opt_parser)

# ── Validate arguments ────────────────────────────────────────────────────────

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: Input file (-i/--input) is required", call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Error: Input file not found:", opt$input), call. = FALSE)
}

if (opt$prevalence < 0 || opt$prevalence > 1) {
  stop("Error: --prevalence must be between 0 and 1", call. = FALSE)
}

valid_dimred <- c("pca", "ica", "rpca")
if (!(opt$dimred %in% valid_dimred)) {
  stop(paste("Error: --dimred must be one of:", paste(valid_dimred, collapse = ", ")),
       call. = FALSE)
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("Error: --cutoff must be strictly between 0 and 1", call. = FALSE)
}

# ── Assign parameters ─────────────────────────────────────────────────────────

input_csv      <- opt$input
output_prefix  <- opt$output
prevalence_thr <- opt$prevalence
min_samp_reads <- opt$minsamp
max_taxa       <- opt$maxtaxa
dim_red_method <- opt$dimred
ncomp          <- opt$ncomp
cutoff         <- opt$cutoff
is_transformed <- opt$transformed
seed           <- opt$seed

set.seed(seed)

# ── Header ────────────────────────────────────────────────────────────────────

message("===================================================")
message("LUPINE Single Time Point Network Inference")
message("===================================================")
message(paste("Input file       :", input_csv))
message(paste("Output prefix    :", output_prefix))
message(paste("Prevalence filter:", prevalence_thr))
message(paste("Min sample reads :", min_samp_reads))
message(paste("Max taxa         :", max_taxa))
message(paste("Dim. reduction   :", dim_red_method, "(ncomp =", ncomp, ")"))
message(paste("p-value cutoff   :", cutoff))
message(paste("Is transformed   :", is_transformed))
message(paste("Random seed      :", seed))
message("")

# ── Load data ─────────────────────────────────────────────────────────────────

message("Loading data...")
data_raw <- read.csv(input_csv, row.names = 1, check.names = FALSE)

message(paste("  Raw dimensions:", nrow(data_raw), "samples x", ncol(data_raw), "taxa"))

data_matrix <- as.matrix(data_raw)
mode(data_matrix) <- "numeric"

# Handle missing values
if (any(is.na(data_matrix))) {
  n_missing <- sum(is.na(data_matrix))
  message(paste("  Warning:", n_missing, "missing values replaced with 0."))
  data_matrix[is.na(data_matrix)] <- 0
}

if (any(data_matrix < 0)) {
  stop("Error: Negative values detected. LUPINE requires non-negative counts.")
}

message("")

# ── Sample filtering ──────────────────────────────────────────────────────────

message("Filtering samples...")
samp_totals     <- rowSums(data_matrix)
keep_samples    <- samp_totals >= min_samp_reads
n_samp_before   <- nrow(data_matrix)
n_samp_after    <- sum(keep_samples)

message(paste("  Min total reads threshold:", min_samp_reads))
message(paste("  Samples before filtering:", n_samp_before))
message(paste("  Samples after  filtering:", n_samp_after))
message(paste("  Samples removed:", n_samp_before - n_samp_after))

if (n_samp_after < 5) {
  stop(paste0(
    "Error: Only ", n_samp_after, " samples remain after filtering. ",
    "LUPINE requires at least 5 samples for reliable partial correlations. ",
    "Lower --minsamp to retain more samples."
  ))
}

data_matrix <- data_matrix[keep_samples, , drop = FALSE]
message("")

# ── Taxa filtering ────────────────────────────────────────────────────────────

message("Filtering taxa...")

# Step 1 – prevalence filter
prevalence      <- colSums(data_matrix > 0) / nrow(data_matrix)
keep_prev       <- prevalence >= prevalence_thr
n_taxa_before   <- ncol(data_matrix)
n_after_prev    <- sum(keep_prev)

message(paste("  Prevalence threshold:", prevalence_thr,
              paste0("(", 100 * prevalence_thr, "% of samples)")))
message(paste("  Taxa before prevalence filter:", n_taxa_before))
message(paste("  Taxa after  prevalence filter:", n_after_prev))

if (n_after_prev == 0) {
  stop("Error: No taxa remain after prevalence filtering. Lower --prevalence.")
}

data_matrix <- data_matrix[, keep_prev, drop = FALSE]

# Step 2 – top-N by mean abundance
if (ncol(data_matrix) > max_taxa) {
  mean_ab   <- colMeans(data_matrix)
  top_idx   <- order(mean_ab, decreasing = TRUE)[seq_len(max_taxa)]
  data_matrix <- data_matrix[, top_idx, drop = FALSE]
  message(paste("  Reduced to top", max_taxa, "taxa by mean abundance"))
}

# Step 3 – remove zero-variance taxa (LUPINE requirement)
taxon_var  <- apply(data_matrix, 2, var)
keep_var   <- taxon_var > 0
n_zero_var <- sum(!keep_var)

if (n_zero_var > 0) {
  message(paste("  Removing", n_zero_var, "zero-variance taxa"))
  data_matrix <- data_matrix[, keep_var, drop = FALSE]
}

n_taxa_final <- ncol(data_matrix)
message(paste("  Final taxa count:", n_taxa_final))

if (n_taxa_final < 2) {
  stop("Error: Fewer than 2 taxa remain after filtering. Adjust filter parameters.")
}

message("")

# ── Summary stats ─────────────────────────────────────────────────────────────

message("Data summary (after filtering):")
message(paste("  Dimensions:", nrow(data_matrix), "samples x", ncol(data_matrix), "taxa"))
message(paste("  Min value:", min(data_matrix)))
message(paste("  Max value:", max(data_matrix)))
message(paste("  Mean value:", round(mean(data_matrix), 2)))
message(paste("  Zero fraction:",
              paste0(round(100 * mean(data_matrix == 0), 1), "%")))
message("")

# ── Build 3D array for LUPINE ─────────────────────────────────────────────────
#
# LUPINE_single() expects a 3D array: samples × taxa × time_points.
# For a single time point we add a singleton third dimension.

message("Preparing data for LUPINE (wrapping 2D matrix into 3D array)...")

data_3d <- array(
  data_matrix,
  dim      = c(nrow(data_matrix), ncol(data_matrix), 1),
  dimnames = list(rownames(data_matrix), colnames(data_matrix), "T1")
)

# Library sizes: row sums per sample × time_point (1-column matrix)
if (!is_transformed) {
  lib_size_mat <- matrix(rowSums(data_matrix), ncol = 1,
                         dimnames = list(rownames(data_matrix), "T1"))
  # Guard against zero library sizes
  if (any(lib_size_mat == 0)) {
    warning("Some library sizes are 0 after filtering; replacing with 1 to avoid log(0).")
    lib_size_mat[lib_size_mat == 0] <- 1
  }
} else {
  lib_size_mat <- NULL
}

message("")

# ── Run LUPINE single time point ──────────────────────────────────────────────

message("Running LUPINE (single time point)...")
message("Parameters:")
message(paste("  - dim. reduction method:", dim_red_method))
message(paste("  - ncomp               :", ncomp))
message(paste("  - p-value cutoff      :", cutoff))
message(paste("  - is.transformed      :", is_transformed))
message("")

lupine_result <- LUPINE_single(
  data           = data_3d,
  day            = 1,
  excluded_taxa  = NULL,
  is.transformed = is_transformed,
  lib_size       = lib_size_mat,
  method         = dim_red_method,
  ncomp          = ncomp
)

message("LUPINE completed successfully!")
message("")

# ── Extract outputs ───────────────────────────────────────────────────────────
#
# lupine_result$Estimate  – partial correlation matrix (taxa × taxa)
# lupine_result$pvalue    – p-value matrix             (taxa × taxa)

asso_mat <- lupine_result$Estimate   # partial correlations
pval_mat <- lupine_result$pvalue     # p-values

# Adjacency: 1 where p-value < cutoff AND not on diagonal
adj_mat              <- (pval_mat < cutoff) * 1
adj_mat[is.na(adj_mat)] <- 0
diag(adj_mat)        <- 0

# ── Network statistics ────────────────────────────────────────────────────────

n_nodes   <- nrow(adj_mat)
n_edges   <- sum(adj_mat, na.rm = TRUE) / 2
max_edges <- n_nodes * (n_nodes - 1) / 2
density   <- if (max_edges > 0) n_edges / max_edges else 0

node_degrees <- rowSums(adj_mat, na.rm = TRUE)

message("Network statistics:")
message(paste("  Nodes (taxa)        :", n_nodes))
message(paste("  Edges               :", n_edges))
message(paste("  Network density     :", round(density, 4)))
message(paste("  Avg degree          :", round(mean(node_degrees), 2)))
message(paste("  Max degree          :", max(node_degrees)))
message(paste("  Isolated nodes      :", sum(node_degrees == 0)))

# Positive / negative edge split
if (!is.null(asso_mat) && n_edges > 0) {
  edge_idx   <- which(adj_mat == 1 & upper.tri(adj_mat))
  asso_vals  <- asso_mat[edge_idx]
  pos_edges  <- sum(asso_vals > 0, na.rm = TRUE)
  message(paste("  Positive edges      :", pos_edges, "/", n_edges,
                paste0("(", round(100 * pos_edges / n_edges, 1), "%)")))
}

message("")

# ── Save results ──────────────────────────────────────────────────────────────

message("Saving results...")

# 1. RData – full LUPINE result object (mirrors net_construct naming used by
#    compare_association_methods.R; stored as 'lupine_result' + helper matrices)
rdata_file <- paste0(output_prefix, "_lupine.RData")
save(lupine_result, asso_mat, pval_mat, adj_mat, data_matrix,
     file = rdata_file)
message(paste("  - Full results saved to:", rdata_file))

# 2. Association (partial correlation) matrix
if (!is.null(asso_mat)) {
  asso_file <- paste0(output_prefix, "_lupine_association_matrix.csv")
  write.csv(asso_mat, file = asso_file)
  message(paste("  - Association matrix saved to:", asso_file))
}

# 3. Adjacency matrix
adj_file <- paste0(output_prefix, "_lupine_adjacency_matrix.csv")
write.csv(adj_mat, file = adj_file)
message(paste("  - Adjacency matrix saved to:", adj_file))

# 4. Edge list
edge_idx_upper <- which(adj_mat == 1 & upper.tri(adj_mat), arr.ind = TRUE)

if (nrow(edge_idx_upper) > 0) {
  edge_list <- data.frame(
    from        = rownames(adj_mat)[edge_idx_upper[, 1]],
    to          = colnames(adj_mat)[edge_idx_upper[, 2]],
    association = if (!is.null(asso_mat)) asso_mat[edge_idx_upper] else NA_real_,
    pvalue      = pval_mat[edge_idx_upper]
  )
  edge_file <- paste0(output_prefix, "_lupine_edge_list.csv")
  write.csv(edge_list, file = edge_file, row.names = FALSE)
  message(paste("  - Edge list saved to:", edge_file))
} else {
  message("  - No significant edges found; edge list not written.")
}

# 5. Node attributes
node_attrs <- data.frame(
  taxon          = colnames(data_matrix),
  degree         = node_degrees,
  mean_abundance = colMeans(data_matrix),
  prevalence     = colSums(data_matrix > 0) / nrow(data_matrix)
)
node_file <- paste0(output_prefix, "_lupine_node_attributes.csv")
write.csv(node_attrs, file = node_file, row.names = FALSE)
message(paste("  - Node attributes saved to:", node_file))

# 6. P-value matrix (extra output useful for downstream inspection)
pval_file <- paste0(output_prefix, "_lupine_pvalue_matrix.csv")
write.csv(pval_mat, file = pval_file)
message(paste("  - P-value matrix saved to:", pval_file))

message("")
message("===================================================")
message("LUPINE single time point analysis completed!")
message("===================================================")
