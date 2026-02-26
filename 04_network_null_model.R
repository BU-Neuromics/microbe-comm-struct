#!/usr/bin/env Rscript
# =============================================================================
# Script: 04_network_null_model.R
#
# Description:
#   Null model validation of the FastSpar co-occurrence network. Builds an
#   igraph network from the filtered edge list, computes observed network
#   metrics (modularity, clustering coefficient, mean path length, degree
#   distribution, cohesion), and compares against 1000 Erdos-Renyi and 1000
#   configuration model null graphs. Computes z-scores and empirical p-values.
#
# Inputs (positional):
#   1. edges_csv  - Filtered edge list CSV from Script 3 (fastspar_edges.csv)
#   2. counts_csv - Counts matrix CSV (rows = taxa, cols = samples); used for
#                   node annotation (library-size independent)
#
# Outputs (saved to same directory as edges_csv):
#   - null_model_comparison.pdf - Two-panel plot: modularity and clustering
#                                 coefficient vs null distributions
#
# Stdout:
#   - Observed network metrics
#   - Z-scores and empirical p-values vs both null models
#   - Module sizes
#   - Network cohesion metrics (Herren & McMahon 2018)
#
# Dependencies:
#   igraph, ggplot2, ggpubr (or patchwork for multi-panel)
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(parallel)
})

# --- Parse arguments ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript 04_network_null_model.R <edges_csv> <counts_csv> [ncores]\n",
      file = stderr())
  quit(status = 1)
}

edges_path  <- args[1]
counts_path <- args[2]
out_dir     <- dirname(edges_path)

# Optional 3rd positional arg: number of cores for null-model simulation loops.
if (length(args) >= 3) {
  n_cores <- max(1L, suppressWarnings(as.integer(args[3])))
  if (is.na(n_cores)) n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
} else {
  n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
}
message("Using ", n_cores, " core(s) for null model simulations.")

# --- Load data ----------------------------------------------------------------
message("Loading edge list from: ", edges_path)
edges_df <- read.csv(edges_path, stringsAsFactors = FALSE)

required_cols <- c("taxon1", "taxon2", "correlation")
missing_req   <- setdiff(required_cols, colnames(edges_df))
if (length(missing_req) > 0) {
  stop("Edge list missing required columns: ",
       paste(missing_req, collapse = ", "), call. = FALSE)
}
message("Loaded ", nrow(edges_df), " edges.")

if (nrow(edges_df) == 0) {
  stop("Edge list is empty. Cannot build network.", call. = FALSE)
}

message("Loading counts matrix from: ", counts_path)
counts_raw <- read.csv(counts_path, row.names = 1, check.names = FALSE)

# --- Build igraph network -----------------------------------------------------
message("Building igraph network...")

graph_obs <- graph_from_data_frame(
  d = data.frame(
    from   = edges_df$taxon1,
    to     = edges_df$taxon2,
    weight = abs(edges_df$correlation),
    sign   = ifelse(edges_df$correlation > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  ),
  directed = FALSE
)

n_nodes <- vcount(graph_obs)
n_edges <- ecount(graph_obs)
message("Network: ", n_nodes, " nodes, ", n_edges, " edges")

# --- Helper: compute network metrics ------------------------------------------
compute_metrics <- function(g) {
  # Modularity via Louvain community detection
  comm       <- cluster_louvain(g)
  mod        <- modularity(comm)

  # Global clustering coefficient (transitivity)
  cc         <- transitivity(g, type = "global")

  # Mean path length on largest connected component
  components <- components(g)
  lcc_nodes  <- which(components$membership ==
                        which.max(components$csize))
  g_lcc      <- induced_subgraph(g, lcc_nodes)
  if (vcount(g_lcc) > 1) {
    mpl <- mean_distance(g_lcc, directed = FALSE)
  } else {
    mpl <- NA_real_
  }

  list(modularity = mod, clustering_coef = cc, mean_path_length = mpl,
       community = comm)
}

# --- Observed metrics ---------------------------------------------------------
message("Computing observed network metrics...")
obs_metrics <- compute_metrics(graph_obs)

obs_mod <- obs_metrics$modularity
obs_cc  <- obs_metrics$clustering_coef
obs_mpl <- obs_metrics$mean_path_length

# Degree distribution
deg_obs      <- degree(graph_obs)
prop_negative <- sum(edges_df$correlation < 0) / nrow(edges_df)

# Module sizes
comm_obs    <- obs_metrics$community
module_sizes <- sort(sizes(comm_obs), decreasing = TRUE)
n_modules    <- length(module_sizes)

cat("\n=== Observed Network Metrics ===\n")
cat("Nodes:", n_nodes, "\n")
cat("Edges:", n_edges, "\n")
cat("Modularity (Louvain):", round(obs_mod, 4), "\n")
cat("Global clustering coefficient:", round(obs_cc, 4), "\n")
cat("Mean path length (LCC):", round(obs_mpl, 4), "\n")
cat("Number of modules:", n_modules, "\n")
cat("Proportion negative edges:", round(prop_negative, 4), "\n")
cat("Module sizes (all):", paste(module_sizes, collapse = ", "), "\n")
if (n_modules >= 3) {
  cat("Three largest modules:", module_sizes[1], module_sizes[2],
      module_sizes[3], "\n")
}
cat("\n")

# --- Cohesion metrics (Herren & McMahon 2018) ---------------------------------
# Positive cohesion: mean of positive correlations (weighted by their values)
# Negative cohesion: mean of absolute values of negative correlations
pos_corrs <- edges_df$correlation[edges_df$correlation > 0]
neg_corrs <- edges_df$correlation[edges_df$correlation < 0]

pos_cohesion <- if (length(pos_corrs) > 0) mean(pos_corrs)      else 0
neg_cohesion <- if (length(neg_corrs) > 0) mean(abs(neg_corrs)) else 0

cat("=== Network Cohesion (Herren & McMahon 2018) ===\n")
cat("Positive cohesion:", round(pos_cohesion, 4), "\n")
cat("Negative cohesion:", round(neg_cohesion, 4), "\n\n")

# --- Null model simulations ---------------------------------------------------
n_null <- 1000
message("Generating ", n_null, " Erdos-Renyi null graphs (", n_cores, " core(s))...")

er_null_results <- mclapply(seq_len(n_null), function(i) {
  set.seed(42L + i)
  g_er <- erdos.renyi.game(n_nodes, n_edges, type = "gnm", directed = FALSE)
  comm_er <- cluster_louvain(g_er)
  list(mod = modularity(comm_er),
       cc  = transitivity(g_er, type = "global"))
}, mc.cores = n_cores)

er_mod_null <- vapply(er_null_results, `[[`, numeric(1), "mod")
er_cc_null  <- vapply(er_null_results, `[[`, numeric(1), "cc")

message("Generating ", n_null, " configuration model null graphs (", n_cores, " core(s))...")
deg_seq   <- deg_obs  # degree sequence from observed network

cm_null_results <- mclapply(seq_len(n_null), function(i) {
  set.seed(42L + n_null + i)
  g_cm <- tryCatch(
    sample_degseq(deg_seq, method = "vl"),
    error = function(e) {
      tryCatch(sample_degseq(deg_seq, method = "simple"),
               error = function(e2) NULL)
    }
  )
  if (is.null(g_cm)) return(list(mod = NA_real_, cc = NA_real_))
  comm_cm <- cluster_louvain(g_cm)
  list(mod = modularity(comm_cm),
       cc  = transitivity(g_cm, type = "global"))
}, mc.cores = n_cores)

cm_mod_null <- vapply(cm_null_results, `[[`, numeric(1), "mod")
cm_cc_null  <- vapply(cm_null_results, `[[`, numeric(1), "cc")

# Remove NAs for statistics
er_mod_null_clean <- er_mod_null[!is.na(er_mod_null)]
er_cc_null_clean  <- er_cc_null[!is.na(er_cc_null)]
cm_mod_null_clean <- cm_mod_null[!is.na(cm_mod_null)]
cm_cc_null_clean  <- cm_cc_null[!is.na(cm_cc_null)]

# --- Z-scores and empirical p-values ------------------------------------------
z_score <- function(obs, null_vec) {
  (obs - mean(null_vec)) / sd(null_vec)
}
emp_pval <- function(obs, null_vec, alternative = "greater") {
  if (alternative == "greater") {
    (sum(null_vec >= obs) + 1) / (length(null_vec) + 1)
  } else {
    (sum(null_vec <= obs) + 1) / (length(null_vec) + 1)
  }
}

cat("=== Null Model Comparison: Modularity ===\n")
z_mod_er  <- z_score(obs_mod, er_mod_null_clean)
p_mod_er  <- emp_pval(obs_mod, er_mod_null_clean, "greater")
z_mod_cm  <- z_score(obs_mod, cm_mod_null_clean)
p_mod_cm  <- emp_pval(obs_mod, cm_mod_null_clean, "greater")

cat("  vs Erdos-Renyi:       z =", round(z_mod_er, 3),
    "  empirical p =", signif(p_mod_er, 3), "\n")
cat("  vs Config model:      z =", round(z_mod_cm, 3),
    "  empirical p =", signif(p_mod_cm, 3), "\n\n")

cat("=== Null Model Comparison: Clustering Coefficient ===\n")
z_cc_er <- z_score(obs_cc, er_cc_null_clean)
p_cc_er <- emp_pval(obs_cc, er_cc_null_clean, "greater")
z_cc_cm <- z_score(obs_cc, cm_cc_null_clean)
p_cc_cm <- emp_pval(obs_cc, cm_cc_null_clean, "greater")

cat("  vs Erdos-Renyi:       z =", round(z_cc_er, 3),
    "  empirical p =", signif(p_cc_er, 3), "\n")
cat("  vs Config model:      z =", round(z_cc_cm, 3),
    "  empirical p =", signif(p_cc_cm, 3), "\n\n")

# --- Plots --------------------------------------------------------------------
message("Generating null model comparison plots...")

# Build data frames for plotting
er_null_df <- data.frame(
  metric = rep(c("Modularity", "Clustering Coefficient"),
               c(length(er_mod_null_clean), length(er_cc_null_clean))),
  value  = c(er_mod_null_clean, er_cc_null_clean),
  model  = "Erdos-Renyi",
  stringsAsFactors = FALSE
)

cm_null_df <- data.frame(
  metric = rep(c("Modularity", "Clustering Coefficient"),
               c(length(cm_mod_null_clean), length(cm_cc_null_clean))),
  value  = c(cm_mod_null_clean, cm_cc_null_clean),
  model  = "Configuration Model",
  stringsAsFactors = FALSE
)

null_df <- rbind(er_null_df, cm_null_df)

obs_line_df <- data.frame(
  metric = c("Modularity", "Clustering Coefficient"),
  obs    = c(obs_mod, obs_cc),
  stringsAsFactors = FALSE
)

# Modularity panel
mod_df <- null_df[null_df$metric == "Modularity", ]
obs_mod_val <- obs_line_df$obs[obs_line_df$metric == "Modularity"]

p_mod <- ggplot(mod_df, aes(x = value, fill = model)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
  geom_vline(xintercept = obs_mod_val, color = "black",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = obs_mod_val, y = Inf,
           label = paste0("Observed\n", round(obs_mod_val, 3)),
           hjust = -0.1, vjust = 1.5, size = 3.5) +
  scale_fill_manual(values = c("Erdos-Renyi" = "steelblue",
                                "Configuration Model" = "tomato")) +
  labs(
    title = "Modularity vs Null Models",
    x     = "Modularity",
    y     = "Count",
    fill  = "Null Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Clustering coefficient panel
cc_df <- null_df[null_df$metric == "Clustering Coefficient", ]
obs_cc_val <- obs_line_df$obs[obs_line_df$metric == "Clustering Coefficient"]

p_cc <- ggplot(cc_df, aes(x = value, fill = model)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
  geom_vline(xintercept = obs_cc_val, color = "black",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = obs_cc_val, y = Inf,
           label = paste0("Observed\n", round(obs_cc_val, 3)),
           hjust = -0.1, vjust = 1.5, size = 3.5) +
  scale_fill_manual(values = c("Erdos-Renyi" = "steelblue",
                                "Configuration Model" = "tomato")) +
  labs(
    title = "Clustering Coefficient vs Null Models",
    x     = "Clustering Coefficient",
    y     = "Count",
    fill  = "Null Model"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Combine using base graphics layout (avoid patchwork dependency)
plot_out <- file.path(out_dir, "null_model_comparison.pdf")
pdf(plot_out, width = 14, height = 6)

# Use grid arrangement via ggplot2 internals
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) {
  gridExtra::grid.arrange(p_mod, p_cc, ncol = 2)
} else {
  # Fallback: print on separate pages
  print(p_mod)
  print(p_cc)
}
dev.off()
message("Null model comparison plot saved to: ", plot_out)

message("Script 04_network_null_model.R complete.")
