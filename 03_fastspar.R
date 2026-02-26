#!/usr/bin/env Rscript
# =============================================================================
# Script: 03_fastspar.R
#
# Description:
#   FastSpar correlation analysis pipeline for microbial co-occurrence.
#   Prepares input for the FastSpar command-line tool, runs correlation
#   estimation with bootstrapping for empirical p-values, parses outputs,
#   applies FDR correction, and produces filtered edge lists and heatmap.
#
# Inputs (positional):
#   1. counts_csv   - Path to counts matrix CSV (rows = taxa, cols = samples)
#   2. metadata_csv - Path to metadata CSV (rows = samples, must have 'group')
#   3. taxonomy_csv - (Optional) Path to taxonomy CSV with columns:
#                     taxon_id, kingdom, phylum, class, order, family, genus
#
# Outputs (saved to same directory as counts_csv):
#   - fastspar_edges.csv     - Filtered edge list (FDR < 0.05, |r| > 0.3)
#   - fastspar_heatmap.pdf   - Correlation heatmap of top 50 taxa by degree
#
# Stdout:
#   - Summary: total significant edges, proportion pos/neg, mean/SD correlations
#   - Cross-kingdom edge breakdown (if taxonomy provided)
#
# Dependencies:
#   pheatmap, ggplot2
#   External: fastspar (must be on PATH)
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
})

# --- Parse arguments ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 03_fastspar.R <counts_csv> <metadata_csv> [taxonomy_csv] [ncores]\n",
      file = stderr())
  quit(status = 1)
}

counts_path   <- args[1]
metadata_path <- args[2]

# Detect optional taxonomy (3rd arg) and ncores (3rd or 4th arg):
# ncores may be passed as args[3] when no taxonomy is provided, or args[4] when
# taxonomy is supplied.
taxonomy_path <- NULL
n_cores       <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

if (length(args) >= 3) {
  # If args[3] looks like an integer it is the ncores value, not a file path.
  arg3_as_int <- suppressWarnings(as.integer(args[3]))
  if (!is.na(arg3_as_int)) {
    n_cores <- max(1L, arg3_as_int)
  } else {
    taxonomy_path <- args[3]
    if (length(args) >= 4) {
      arg4_as_int <- suppressWarnings(as.integer(args[4]))
      if (!is.na(arg4_as_int)) {
        n_cores <- max(1L, arg4_as_int)
      }
    }
  }
}

message("Using ", n_cores, " thread(s) for FastSpar.")
out_dir       <- dirname(counts_path)

# --- Check FastSpar is available ----------------------------------------------
fastspar_bin <- Sys.which("fastspar")
if (fastspar_bin == "") {
  stop("FastSpar executable not found on PATH. Please install FastSpar.",
       call. = FALSE)
}
message("FastSpar binary found at: ", fastspar_bin)

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

# --- Load taxonomy (optional) -------------------------------------------------
taxonomy <- NULL
if (!is.null(taxonomy_path)) {
  message("Loading taxonomy from: ", taxonomy_path)
  taxonomy <- read.csv(taxonomy_path, stringsAsFactors = FALSE)
  required_cols <- c("taxon_id", "kingdom", "phylum", "class",
                     "order", "family", "genus")
  missing_cols <- setdiff(required_cols, colnames(taxonomy))
  if (length(missing_cols) > 0) {
    warning("Taxonomy file missing columns: ",
            paste(missing_cols, collapse = ", "),
            ". Cross-kingdom annotation will be partial.")
  }
}

# --- Prepare FastSpar input ---------------------------------------------------
# FastSpar format: TSV, samples as rows, taxa as columns
# First column header is "#OTU ID" (actually the sample ID column label)
# But FastSpar expects: taxa as rows, samples as columns with first row being
# header "#OTU ID\tsample1\tsample2..."
# Verify: FastSpar input is OTU table where rows=OTUs, cols=samples
# Header: #OTU ID <tab> sample1 <tab> sample2 ...

message("Writing FastSpar input TSV...")
tmp_dir     <- tempdir()
fs_tmp      <- file.path(tmp_dir, "fastspar_run")
dir.create(fs_tmp, showWarnings = FALSE, recursive = TRUE)

fs_input    <- file.path(fs_tmp, "counts.tsv")
fs_corr_out <- file.path(fs_tmp, "correlations.tsv")
fs_cov_out  <- file.path(fs_tmp, "covariance.tsv")
bs_dir      <- file.path(fs_tmp, "bootstraps")
bs_corr_dir <- file.path(fs_tmp, "bootstrap_correlations")
pval_out    <- file.path(fs_tmp, "pvalues.tsv")

dir.create(bs_dir,      showWarnings = FALSE)
dir.create(bs_corr_dir, showWarnings = FALSE)

# Build the TSV: first column = taxon IDs, rest = sample columns
# Round to integers required by some FastSpar versions; since our data is
# fractional we round to avoid issues and warn user.
counts_integer <- round(counts_raw)
any_fractional <- any(counts_raw != counts_integer)
if (any_fractional) {
  message("WARNING: Counts contain non-integer values. FastSpar requires ",
          "integer counts. Rounding to nearest integer. Values near 0 may ",
          "become 0.")
}

header_line <- paste(c("#OTU ID", colnames(counts_integer)), collapse = "\t")
row_lines   <- apply(counts_integer, 1, function(row) {
  paste(c(rownames(counts_integer)[which(rownames(counts_integer) ==
                                          names(row)[1])],
          as.integer(row)), collapse = "\t")
})
# Fix: write properly using write.table
counts_df_out <- as.data.frame(counts_integer)
counts_df_out <- cbind("#OTU ID" = rownames(counts_df_out), counts_df_out)

write.table(counts_df_out, file = fs_input, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# FastSpar uses tab-separated; fix the header to use "#OTU ID" not "X.OTU.ID"
lines_in <- readLines(fs_input)
lines_in[1] <- sub("^#OTU\\.ID", "#OTU ID", lines_in[1])
writeLines(lines_in, fs_input)

message("FastSpar input written to: ", fs_input)

# --- Run FastSpar main analysis -----------------------------------------------
n_iter      <- 500
n_excl_iter <- 10

message("Running FastSpar (", n_iter, " iterations, ",
        n_excl_iter, " exclusion iterations)...")

fs_cmd <- paste(
  shQuote(fastspar_bin),
  "--otu_table", shQuote(fs_input),
  "--correlation", shQuote(fs_corr_out),
  "--covariance", shQuote(fs_cov_out),
  "--iterations", n_iter,
  "--exclude_iterations", n_excl_iter,
  "--threads", n_cores,
  "--yes"
)
message("Command: ", fs_cmd)
fs_exit <- system(fs_cmd)
if (fs_exit != 0) {
  stop("FastSpar failed with exit code: ", fs_exit, call. = FALSE)
}
message("FastSpar main analysis complete.")

# --- Run FastSpar bootstrapping -----------------------------------------------
n_bootstraps <- 1000
message("Generating ", n_bootstraps, " bootstrap samples...")

bs_cmd <- paste(
  "fastspar_bootstrap",
  "--otu_table", shQuote(fs_input),
  "--number", n_bootstraps,
  "--prefix", shQuote(file.path(bs_dir, "bootstrap"))
)
# Try fastspar_bootstrap; fall back to fastspar --resample if not available
bs_bin <- Sys.which("fastspar_bootstrap")
if (bs_bin == "") {
  message("fastspar_bootstrap not found; trying fastspar --resample_data...")
  bs_cmd <- paste(
    shQuote(fastspar_bin),
    "--otu_table", shQuote(fs_input),
    "--resample_data", n_bootstraps,
    "--resample_prefix", shQuote(file.path(bs_dir, "bootstrap"))
  )
} else {
  bs_cmd <- paste(
    shQuote(bs_bin),
    "--otu_table", shQuote(fs_input),
    "--number", n_bootstraps,
    "--prefix", shQuote(file.path(bs_dir, "bootstrap"))
  )
}
message("Bootstrap command: ", bs_cmd)
bs_exit <- system(bs_cmd)
if (bs_exit != 0) {
  stop("FastSpar bootstrap generation failed with exit code: ", bs_exit,
       call. = FALSE)
}

message("Computing correlations on each bootstrap sample...")
bs_files <- list.files(bs_dir, pattern = "\\.tsv$", full.names = TRUE)
message("  Found ", length(bs_files), " bootstrap samples")

# Run bootstrap correlations in parallel: n_cores single-threaded fastspar jobs
# run simultaneously rather than one n_cores-threaded job sequentially.
# This is substantially faster, especially under emulation environments.
message("  Dispatching ", n_cores, " parallel workers (--threads 1 each)...")
bs_results <- parallel::mclapply(bs_files, function(bsf) {
  bsf_base  <- sub("\\.tsv$", "", basename(bsf))
  bs_corr_f <- file.path(bs_corr_dir, paste0(bsf_base, "_corr.tsv"))
  bs_cov_f  <- file.path(bs_corr_dir, paste0(bsf_base, "_cov.tsv"))
  bs_run_cmd <- paste(
    shQuote(fastspar_bin),
    "--otu_table", shQuote(bsf),
    "--correlation", shQuote(bs_corr_f),
    "--covariance", shQuote(bs_cov_f),
    "--iterations", 200,
    "--exclude_iterations", 5,
    "--threads", 1,
    "--yes"
  )
  run_exit <- system(bs_run_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  list(bsf = bsf, exit = run_exit)
}, mc.cores = n_cores)

failed_bs <- Filter(function(r) r$exit != 0, bs_results)
if (length(failed_bs) > 0) {
  warning(length(failed_bs), " bootstrap run(s) failed. Tmp dir: ", fs_tmp)
}
message("Bootstrap correlation computation complete.")

message("Computing empirical p-values from bootstrap distribution...")
pval_bin <- Sys.which("fastspar_pvalues")
if (pval_bin == "") pval_bin <- Sys.which("fastspar_p_values")
if (pval_bin == "") {
  message("  fastspar_pvalues not found; trying fastspar --p_values subcommand...")
  pval_cmd <- paste(
    shQuote(fastspar_bin),
    "--otu_table", shQuote(fs_input),
    "--correlation", shQuote(fs_corr_out),
    "--prefix", shQuote(file.path(bs_corr_dir, "bootstrap")),
    "--permutations", n_bootstraps,
    "--outfile", shQuote(pval_out)
  )
} else {
  pval_cmd <- paste(
    shQuote(pval_bin),
    "--otu_table", shQuote(fs_input),
    "--correlation", shQuote(fs_corr_out),
    "--prefix", shQuote(file.path(bs_corr_dir, "bootstrap")),
    "--permutations", n_bootstraps,
    "--outfile", shQuote(pval_out)
  )
}
message("P-value command: ", pval_cmd)
pval_exit <- system(pval_cmd)
if (pval_exit != 0) {
  stop("FastSpar p-value computation failed with exit code: ", pval_exit,
       call. = FALSE)
}
message("P-value computation complete.")

# --- Parse outputs ------------------------------------------------------------
message("Parsing FastSpar outputs...")

parse_fastspar_matrix <- function(path) {
  mat <- read.table(path, header = TRUE, sep = "\t", row.names = 1,
                    check.names = FALSE, comment.char = "")
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat
}

corr_mat <- parse_fastspar_matrix(fs_corr_out)
pval_mat <- parse_fastspar_matrix(pval_out)

n_taxa_corr <- nrow(corr_mat)
message("Parsed correlation matrix: ", n_taxa_corr, " x ", ncol(corr_mat))

# Symmetrize (FastSpar output should be symmetric but verify)
corr_mat <- (corr_mat + t(corr_mat)) / 2
pval_mat <- (pval_mat + t(pval_mat)) / 2

# --- FDR correction -----------------------------------------------------------
message("Applying BH FDR correction...")
# Extract upper triangle (excluding diagonal)
upper_idx  <- which(upper.tri(pval_mat))
pvals_vec  <- pval_mat[upper_idx]
fdr_vec    <- p.adjust(pvals_vec, method = "BH")

taxa_names <- rownames(corr_mat)
idx_pairs  <- which(upper.tri(pval_mat), arr.ind = TRUE)

edge_df <- data.frame(
  taxon1      = taxa_names[idx_pairs[, 1]],
  taxon2      = taxa_names[idx_pairs[, 2]],
  correlation = corr_mat[upper_idx],
  pvalue      = pvals_vec,
  fdr         = fdr_vec,
  stringsAsFactors = FALSE
)

# --- Filter edges -------------------------------------------------------------
sig_mask  <- edge_df$fdr < 0.05 & abs(edge_df$correlation) > 0.3
sig_edges <- edge_df[sig_mask, ]
sig_edges <- sig_edges[order(sig_edges$fdr), ]

message("Significant edges (FDR < 0.05, |r| > 0.3): ", nrow(sig_edges))

edges_out <- file.path(out_dir, "fastspar_edges.csv")
write.csv(sig_edges, edges_out, row.names = FALSE)
message("Edge list saved to: ", edges_out)

# --- Summary statistics -------------------------------------------------------
cat("\n=== FastSpar Network Summary ===\n")
cat("Total significant edges:", nrow(sig_edges), "\n")

if (nrow(sig_edges) > 0) {
  n_pos  <- sum(sig_edges$correlation > 0)
  n_neg  <- sum(sig_edges$correlation < 0)
  cat("Positive edges:", n_pos, "(", round(100 * n_pos / nrow(sig_edges), 1), "%)\n")
  cat("Negative edges:", n_neg, "(", round(100 * n_neg / nrow(sig_edges), 1), "%)\n")
  cat("Mean |correlation|:", round(mean(abs(sig_edges$correlation)), 4), "\n")
  cat("SD   |correlation|:", round(sd(abs(sig_edges$correlation)), 4), "\n\n")
}

# --- Cross-kingdom breakdown (if taxonomy provided) ---------------------------
if (!is.null(taxonomy) && nrow(sig_edges) > 0 &&
    all(c("taxon_id", "kingdom") %in% colnames(taxonomy))) {
  taxon_kingdom <- setNames(taxonomy$kingdom, taxonomy$taxon_id)
  sig_edges$kingdom1 <- taxon_kingdom[sig_edges$taxon1]
  sig_edges$kingdom2 <- taxon_kingdom[sig_edges$taxon2]

  # Create pair label (alphabetically sorted)
  kingdom_pair <- apply(sig_edges[, c("kingdom1", "kingdom2")], 1, function(r) {
    k <- sort(r)
    paste(k, collapse = " <-> ")
  })
  kingdom_table <- table(kingdom_pair)
  cat("=== Cross-kingdom Edge Breakdown ===\n")
  print(as.data.frame(kingdom_table))
  cat("\n")
}

# --- Heatmap: top 50 taxa by degree -------------------------------------------
message("Generating correlation heatmap...")
if (nrow(sig_edges) > 0) {
  # Compute degree for each taxon
  all_taxa_sig <- c(sig_edges$taxon1, sig_edges$taxon2)
  degree_table <- sort(table(all_taxa_sig), decreasing = TRUE)
  top50_taxa   <- names(degree_table)[seq_len(min(50, length(degree_table)))]

  # Subset correlation matrix
  top50_taxa_avail <- intersect(top50_taxa, rownames(corr_mat))
  corr_sub <- corr_mat[top50_taxa_avail, top50_taxa_avail]

  # Annotation columns (kingdom if taxonomy available)
  annotation_row <- NULL
  ann_colors     <- NULL
  if (!is.null(taxonomy) && "kingdom" %in% colnames(taxonomy)) {
    taxon_kingdom_df <- data.frame(
      Kingdom = taxon_kingdom[top50_taxa_avail],
      row.names = top50_taxa_avail,
      stringsAsFactors = FALSE
    )
    annotation_row <- taxon_kingdom_df

    kingdoms <- unique(na.omit(taxon_kingdom_df$Kingdom))
    pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
             "#A65628", "#F781BF", "#999999")
    kingdom_colors <- setNames(pal[seq_along(kingdoms)], kingdoms)
    ann_colors <- list(Kingdom = kingdom_colors)
  }

  heatmap_out <- file.path(out_dir, "fastspar_heatmap.pdf")
  pdf(heatmap_out, width = 14, height = 12)
  pheatmap(
    corr_sub,
    clustering_method   = "complete",
    color               = colorRampPalette(c("steelblue", "white", "tomato"))(100),
    breaks              = seq(-1, 1, length.out = 101),
    annotation_row      = annotation_row,
    annotation_col      = annotation_row,
    annotation_colors   = ann_colors,
    fontsize_row        = 6,
    fontsize_col        = 6,
    main                = "FastSpar Correlations — Top 50 Taxa by Degree",
    show_rownames       = TRUE,
    show_colnames       = TRUE,
    border_color        = NA
  )
  dev.off()
  message("Heatmap saved to: ", heatmap_out)
} else {
  message("No significant edges found; skipping heatmap.")
}

message("Script 03_fastspar.R complete.")
