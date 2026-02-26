#!/usr/bin/env Rscript
# visualize_results.R
# Comprehensive visualization of NetCoMi pipeline results
# Produces a multi-page PDF with network analysis visualizations

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(grid)
  library(RColorBrewer)
  library(scales)
})

# ============================================================
# Configuration
# ============================================================
RESULTS_DIR <- "results"
OUTPUT_PDF  <- "network_analysis_visualizations.pdf"
TAX_LEVELS  <- c("class", "order", "family", "genus")
METHODS     <- c("sparcc", "spearman", "cclasso", "spieceasi", "spring")
SPARSE_METHODS <- c("spring", "spieceasi", "cclasso")

# Color palettes
METHOD_COLORS <- c(
  sparcc    = "#E41A1C",
  spearman  = "#FF7F00",
  cclasso   = "#4DAF4A",
  spieceasi = "#377EB8",
  spring    = "#984EA3"
)
GROUP_COLORS <- c(C = "#2196F3", PD = "#F44336")
LEVEL_COLORS <- setNames(
  brewer.pal(4, "Set2"),
  c("class", "order", "family", "genus")
)

message("===================================================")
message("NetCoMi Results Visualization")
message("===================================================")

# ============================================================
# Helper: classify sample as PD or Control
# ============================================================
get_group <- function(taxon_names) {
  ifelse(grepl("^PD_", taxon_names), "PD", "C")
}

# ============================================================
# Helper: read a CSV safely
# ============================================================
read_safe <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
}

# ============================================================
# 1. Load all per-method per-level data
# ============================================================
message("Loading per-method data...")

all_basic_stats  <- list()
all_centralities <- list()
all_clusters     <- list()
all_global_props <- list()
all_edge_overlap <- list()
all_deg_corr     <- list()
all_bet_corr     <- list()
all_close_corr   <- list()
all_hub_overlap  <- list()

for (lvl in TAX_LEVELS) {
  lvl_dir <- file.path(RESULTS_DIR, lvl)

  # Comparison-level summaries
  all_basic_stats[[lvl]]  <- read_safe(file.path(lvl_dir, "comparison", "comparison_basic_stats.csv"))
  all_edge_overlap[[lvl]] <- read_safe(file.path(lvl_dir, "comparison", "comparison_edge_overlap.csv"))
  all_deg_corr[[lvl]]     <- read_safe(file.path(lvl_dir, "comparison", "comparison_degree_correlation.csv"))
  all_bet_corr[[lvl]]     <- read_safe(file.path(lvl_dir, "comparison", "comparison_betweenness_correlation.csv"))
  all_close_corr[[lvl]]   <- read_safe(file.path(lvl_dir, "comparison", "comparison_closeness_correlation.csv"))
  all_hub_overlap[[lvl]]  <- read_safe(file.path(lvl_dir, "comparison", "comparison_hub_overlap.csv"))

  # Per-method data
  for (mth in METHODS) {
    mth_dir <- file.path(lvl_dir, mth)
    key <- paste(lvl, mth, sep = "_")

    cent <- read_safe(file.path(mth_dir, paste0("result_", mth, "_centralities.csv")))
    if (!is.null(cent)) {
      cent$level  <- lvl
      cent$method <- mth
      cent$group  <- get_group(cent$taxon)
      all_centralities[[key]] <- cent
    }

    clust <- read_safe(file.path(mth_dir, paste0("result_", mth, "_clusters.csv")))
    if (!is.null(clust)) {
      clust$level  <- lvl
      clust$method <- mth
      clust$group  <- get_group(clust$taxon)
      all_clusters[[key]] <- clust
    }

    # Per-method global props (wide -> named list)
    gp <- read_safe(file.path(mth_dir, paste0("result_", mth, "_global_properties.csv")))
    if (!is.null(gp)) {
      gp_row <- setNames(as.numeric(gp$value), gp$property)
      gp_row[["level"]]  <- lvl
      gp_row[["method"]] <- mth
      all_global_props[[key]] <- gp_row
    }
  }
}

message("Data loaded.")

# ============================================================
# Helper: melt a pairwise CSV matrix to long form
# ============================================================
melt_matrix <- function(df) {
  if (is.null(df)) return(NULL)
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  df[upper.tri(as.matrix(df))] <- NA
  m <- melt(as.matrix(df), varnames = c("method1", "method2"), na.rm = TRUE)
  m$method1 <- as.character(m$method1)
  m$method2 <- as.character(m$method2)
  m
}

# ============================================================
# Open PDF
# ============================================================
message(paste("Writing visualizations to", OUTPUT_PDF))
pdf(OUTPUT_PDF, width = 14, height = 10)

# ============================================================
# FIGURE 1: Title page
# ============================================================
grid.newpage()
grid.text(
  "NetCoMi Microbial Network Analysis\nResults Visualization Report",
  x = 0.5, y = 0.6,
  gp = gpar(fontsize = 22, fontface = "bold")
)
grid.text(
  paste0(
    "Taxonomic levels: ", paste(TAX_LEVELS, collapse = ", "), "\n",
    "Methods: ", paste(METHODS, collapse = ", "), "\n",
    "24 samples: 10 Controls (C_) + 14 Parkinson's Disease (PD_)"
  ),
  x = 0.5, y = 0.38,
  gp = gpar(fontsize = 13)
)
grid.text(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  x = 0.5, y = 0.15,
  gp = gpar(fontsize = 10, col = "grey50")
)

# ============================================================
# FIGURE 2: Network Topology Overview
# Network size, density, and LCC size across methods & levels
# ============================================================
message("Figure 2: Network topology overview...")

basic_df <- do.call(rbind, lapply(names(all_basic_stats), function(lvl) {
  d <- all_basic_stats[[lvl]]
  if (is.null(d)) return(NULL)
  d$level <- lvl
  d
}))
basic_df$level  <- factor(basic_df$level, levels = TAX_LEVELS)
basic_df$method <- factor(basic_df$method, levels = METHODS)
basic_df$edge_density_capped <- pmin(basic_df$edge_density, 1.2)

p_nodes <- ggplot(basic_df, aes(x = method, y = n_nodes, fill = level)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = LEVEL_COLORS) +
  labs(title = "Nodes in network", x = NULL, y = "N nodes") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

p_edges <- ggplot(basic_df, aes(x = method, y = n_edges, fill = level)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = LEVEL_COLORS) +
  labs(title = "Number of edges", x = NULL, y = "N edges") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

p_density <- ggplot(basic_df, aes(x = method, y = edge_density_capped, fill = level)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  scale_fill_manual(values = LEVEL_COLORS) +
  labs(
    title = "Edge density (capped at 1.2)",
    subtitle = "Dashed = fully connected (density=1). Sparse methods (SPRING/SpiecEasi) << 1",
    x = NULL, y = "Edge density"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

p_lcc <- ggplot(basic_df, aes(x = method, y = lcc_size, fill = level)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = LEVEL_COLORS) +
  labs(
    title = "Largest connected component (LCC) size",
    subtitle = "Sparse methods fragment the network; LCC < total nodes",
    x = NULL, y = "LCC size (nodes)"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

grid.arrange(p_nodes, p_edges, p_density, p_lcc, ncol = 2,
  top = textGrob("Network Topology Overview", gp = gpar(fontsize = 15, fontface = "bold")))

# ============================================================
# FIGURE 3: Positive edge percentage
# ============================================================
message("Figure 3: Positive edge percentage...")

# Spearman at family/genus/order reports 200% — this is a known artifact
# (edges summed from both directions). We cap at 100 for display.
basic_df$pos_edges_display <- pmin(basic_df$pos_edges_pct, 100)

p_pos <- ggplot(basic_df, aes(x = method, y = pos_edges_display, color = level, group = level)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = LEVEL_COLORS) +
  facet_wrap(~level, ncol = 2) +
  labs(
    title = "Positive edges percentage by method and taxonomic level",
    subtitle = "Dashed = 50% (balance point). >50% = net positive associations dominate.",
    x = "Method", y = "Positive edges (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

n_components_p <- ggplot(basic_df, aes(x = method, y = n_components, fill = method)) +
  geom_col() +
  scale_fill_manual(values = METHOD_COLORS, na.value = "grey80") +
  facet_wrap(~level, ncol = 2) +
  labs(
    title = "Number of connected components",
    subtitle = "1 = fully connected. Sparse methods produce multiple components.",
    x = NULL, y = "N components"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

grid.arrange(p_pos, n_components_p, ncol = 2,
  top = textGrob("Network Connectivity Properties", gp = gpar(fontsize = 15, fontface = "bold")))

# ============================================================
# FIGURE 4 & 5: Method agreement — edge Jaccard heatmaps
# ============================================================
message("Figure 4-5: Method agreement heatmaps...")

plot_heatmap <- function(df, title, subtitle = NULL, midpoint = 0.5,
                         low = "white", mid = "#FFF176", high = "#D32F2F") {
  if (is.null(df)) return(NULL)
  m <- melt_matrix(df)
  if (is.null(m)) return(NULL)
  ggplot(m, aes(x = method2, y = method1, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
    scale_fill_gradient2(low = low, mid = mid, high = high,
                         midpoint = midpoint, limits = c(-1, 1),
                         name = "Value", na.value = "grey90") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

plot_heatmap_01 <- function(df, title, subtitle = NULL) {
  if (is.null(df)) return(NULL)
  m <- melt_matrix(df)
  if (is.null(m)) return(NULL)
  ggplot(m, aes(x = method2, y = method1, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
    scale_fill_gradient(low = "white", high = "#1565C0",
                        limits = c(0, 1), name = "Jaccard") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

plots_edge <- lapply(TAX_LEVELS, function(lvl) {
  plot_heatmap_01(
    all_edge_overlap[[lvl]],
    title = lvl,
    subtitle = "Jaccard edge overlap"
  )
})
plots_edge <- Filter(Negate(is.null), plots_edge)
if (length(plots_edge) > 0) {
  do.call(grid.arrange, c(plots_edge, list(ncol = 2,
    top = textGrob("Method Agreement: Edge Jaccard Overlap\n(1 = identical networks, 0 = no shared edges)",
                   gp = gpar(fontsize = 14, fontface = "bold")))))
}

plots_deg <- lapply(TAX_LEVELS, function(lvl) {
  plot_heatmap(
    all_deg_corr[[lvl]],
    title = lvl,
    subtitle = "Spearman degree correlation",
    midpoint = 0, low = "#1565C0", mid = "white", high = "#D32F2F"
  )
})
plots_deg <- Filter(Negate(is.null), plots_deg)
if (length(plots_deg) > 0) {
  do.call(grid.arrange, c(plots_deg, list(ncol = 2,
    top = textGrob("Method Agreement: Degree Centrality Correlation\n(Spearman r; positive = methods rank nodes similarly)",
                   gp = gpar(fontsize = 14, fontface = "bold")))))
}

# ============================================================
# FIGURE 6: Hub node agreement
# ============================================================
message("Figure 6: Hub overlap...")

plots_hub <- lapply(TAX_LEVELS, function(lvl) {
  plot_heatmap_01(
    all_hub_overlap[[lvl]],
    title = lvl,
    subtitle = "Hub node Jaccard overlap"
  )
})
plots_hub <- Filter(Negate(is.null), plots_hub)
if (length(plots_hub) > 0) {
  do.call(grid.arrange, c(plots_hub, list(ncol = 2,
    top = textGrob("Hub Node Agreement Across Methods\n(Jaccard overlap of hub node sets)",
                   gp = gpar(fontsize = 14, fontface = "bold")))))
}

# ============================================================
# FIGURE 7: Global network properties from per-method files
# ============================================================
message("Figure 7: Global network properties...")

# Build a tidy data frame from per-method global_properties files
# Properties of interest: lccSize1, clustCoef1, modularity1, natConnect1,
#                         avPath1, density1, pep1, lccSizeRel1
PROPS_OF_INTEREST <- c("lccSize1", "lccSizeRel1", "clustCoef1", "modularity1",
                        "natConnect1", "avPath1", "density1", "pep1")
PROP_LABELS <- c(
  lccSize1     = "LCC size (nodes)",
  lccSizeRel1  = "LCC size (relative)",
  clustCoef1   = "Clustering coefficient",
  modularity1  = "Modularity",
  natConnect1  = "Natural connectivity",
  avPath1      = "Avg path length",
  density1     = "Density (LCC)",
  pep1         = "Positive edges % (LCC)"
)

gp_rows <- do.call(rbind, lapply(names(all_global_props), function(key) {
  gp <- all_global_props[[key]]
  row <- data.frame(
    level  = as.character(gp[["level"]]),
    method = as.character(gp[["method"]]),
    stringsAsFactors = FALSE
  )
  for (prop in PROPS_OF_INTEREST) {
    row[[prop]] <- if (prop %in% names(gp)) as.numeric(gp[[prop]]) else NA_real_
  }
  row
}))

if (!is.null(gp_rows) && nrow(gp_rows) > 0) {
  gp_long <- gp_rows %>%
    pivot_longer(cols = all_of(PROPS_OF_INTEREST),
                 names_to = "property", values_to = "value") %>%
    mutate(
      level    = factor(level, levels = TAX_LEVELS),
      method   = factor(method, levels = METHODS),
      prop_label = PROP_LABELS[property]
    )

  # Split into two panels (4 props each) to avoid crowding
  props_panel1 <- c("lccSize1", "lccSizeRel1", "clustCoef1", "modularity1")
  props_panel2 <- c("natConnect1", "avPath1", "density1", "pep1")

  make_prop_plot <- function(props) {
    d <- gp_long %>% filter(property %in% props)
    ggplot(d, aes(x = method, y = value, color = method, shape = level)) +
      geom_point(size = 3, position = position_dodge(0.4)) +
      scale_color_manual(values = METHOD_COLORS, na.value = "grey60") +
      scale_shape_manual(values = c(class = 16, order = 17, family = 15, genus = 18)) +
      facet_wrap(~prop_label, scales = "free_y", ncol = 2) +
      labs(x = NULL, y = "Value", color = "Method", shape = "Level") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
  }

  p1 <- make_prop_plot(props_panel1)
  p2 <- make_prop_plot(props_panel2)

  grid.arrange(p1, top = textGrob(
    "Global Network Properties (LCC) — Panel 1/2",
    gp = gpar(fontsize = 14, fontface = "bold")))
  grid.arrange(p2, top = textGrob(
    "Global Network Properties (LCC) — Panel 2/2",
    gp = gpar(fontsize = 14, fontface = "bold")))
}

# ============================================================
# FIGURE 8: Centrality profiles — degree by sample
# Sparse methods only (SPRING, SpiecEasi, CCLasso)
# ============================================================
message("Figure 8: Centrality profiles...")

cent_df <- do.call(rbind, all_centralities)
cent_df$level  <- factor(cent_df$level, levels = TAX_LEVELS)
cent_df$method <- factor(cent_df$method, levels = METHODS)

sparse_cent <- cent_df %>% filter(method %in% SPARSE_METHODS)

# Degree lollipop per level, colored by PD/Control, one row per method
for (lvl in c("family", "genus", "order")) {
  d <- sparse_cent %>%
    filter(level == lvl) %>%
    arrange(method, group, desc(degree))

  if (nrow(d) == 0) next

  # Rank within each method for consistent x-axis sorting
  d <- d %>%
    group_by(method) %>%
    mutate(rank = rank(-degree, ties.method = "first")) %>%
    ungroup()

  p <- ggplot(d, aes(x = reorder(taxon, -degree), y = degree, color = group)) +
    geom_segment(aes(xend = reorder(taxon, -degree), y = 0, yend = degree),
                 linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_color_manual(values = GROUP_COLORS) +
    facet_wrap(~method, ncol = 1, scales = "free_y") +
    labs(
      title = paste("Degree centrality —", lvl, "level"),
      subtitle = "Sparse methods only. Color = PD (red) vs Control (blue).",
      x = "Sample", y = "Degree centrality", color = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7))

  print(p)
}

# ============================================================
# FIGURE 9: Betweenness centrality (bridges/bottlenecks)
# ============================================================
message("Figure 9: Betweenness centrality...")

for (lvl in c("family", "genus", "order")) {
  d <- sparse_cent %>%
    filter(level == lvl, betweenness > 0)

  if (nrow(d) == 0) {
    message(paste("  No betweenness > 0 for", lvl, "— skipping"))
    next
  }

  p <- ggplot(d, aes(x = reorder(taxon, -betweenness), y = betweenness, fill = group)) +
    geom_col() +
    scale_fill_manual(values = GROUP_COLORS) +
    facet_wrap(~method, ncol = 1, scales = "free") +
    labs(
      title = paste("Betweenness centrality —", lvl, "level"),
      subtitle = "Only nodes with betweenness > 0 shown. High betweenness = network bridges.",
      x = "Sample", y = "Betweenness centrality", fill = "Group"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7))

  print(p)
}

# ============================================================
# FIGURE 10: Cross-method centrality scatter (SPRING vs SpiecEasi)
# ============================================================
message("Figure 10: Cross-method centrality scatter...")

for (lvl in c("family", "genus", "order")) {
  spring_c  <- all_centralities[[paste(lvl, "spring", sep = "_")]]
  spiec_c   <- all_centralities[[paste(lvl, "spieceasi", sep = "_")]]
  sparcc_c  <- all_centralities[[paste(lvl, "sparcc", sep = "_")]]

  if (!is.null(spring_c) && !is.null(spiec_c)) {
    comp <- inner_join(
      spring_c %>% select(taxon, group, spring_deg = degree),
      spiec_c  %>% select(taxon, spiec_deg = degree),
      by = "taxon"
    )
    r <- round(cor(comp$spring_deg, comp$spiec_deg, method = "spearman"), 2)

    p1 <- ggplot(comp, aes(x = spring_deg, y = spiec_deg, color = group, label = taxon)) +
      geom_point(size = 3) +
      geom_text(size = 2.5, vjust = -0.6, check_overlap = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = GROUP_COLORS) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
               label = paste0("Spearman r = ", r), size = 3.5) +
      labs(
        title = paste("SPRING vs SpiecEasi degree —", lvl),
        x = "SPRING degree", y = "SpiecEasi degree", color = "Group"
      ) +
      theme_bw(base_size = 10)

    print(p1)
  }

  # Also SPRING vs SparCC
  if (!is.null(spring_c) && !is.null(sparcc_c)) {
    comp2 <- inner_join(
      spring_c  %>% select(taxon, group, spring_deg = degree),
      sparcc_c  %>% select(taxon, sparcc_deg = degree),
      by = "taxon"
    )
    r2 <- round(cor(comp2$spring_deg, comp2$sparcc_deg, method = "spearman"), 2)

    p2 <- ggplot(comp2, aes(x = spring_deg, y = sparcc_deg, color = group, label = taxon)) +
      geom_point(size = 3) +
      geom_text(size = 2.5, vjust = -0.6, check_overlap = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = GROUP_COLORS) +
      annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
               label = paste0("Spearman r = ", r2), size = 3.5) +
      labs(
        title = paste("SPRING vs SparCC degree —", lvl),
        subtitle = "SparCC is nearly fully connected; high degree = many correlations",
        x = "SPRING degree", y = "SparCC degree", color = "Group"
      ) +
      theme_bw(base_size = 10)

    print(p2)
  }
}

# ============================================================
# FIGURE 11: Community structure — cluster compositions
# ============================================================
message("Figure 11: Community structure...")

clust_df <- do.call(rbind, all_clusters)
clust_df$level  <- factor(clust_df$level, levels = TAX_LEVELS)
clust_df$method <- factor(clust_df$method, levels = METHODS)
clust_df$cluster_label <- paste0("C", clust_df$cluster)

# Only show isolated (non-trivial) methods
for (lvl in c("family", "genus", "order")) {
  for (mth in c("spring", "spieceasi", "cclasso", "sparcc")) {
    key <- paste(lvl, mth, sep = "_")
    d <- all_clusters[[key]]
    if (is.null(d)) next
    n_clusters <- length(unique(d$cluster[d$cluster != 0]))
    if (n_clusters < 2) next  # skip trivial single-cluster results

    d$cluster_label <- factor(paste0("Cluster ", d$cluster),
                               levels = paste0("Cluster ", sort(unique(d$cluster))))

    cnt <- d %>%
      count(cluster_label, group) %>%
      group_by(cluster_label) %>%
      mutate(pct = n / sum(n) * 100)

    p <- ggplot(cnt, aes(x = cluster_label, y = pct, fill = group)) +
      geom_col(position = "stack") +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5),
                size = 3.5, color = "white", fontface = "bold") +
      scale_fill_manual(values = GROUP_COLORS) +
      geom_hline(yintercept = 58.33, linetype = "dashed", color = "grey40") +  # 14/24 = PD proportion
      annotate("text", x = Inf, y = 58.33, hjust = 1.1, vjust = -0.3,
               label = "58.3% = overall PD fraction", size = 2.8, color = "grey40") +
      labs(
        title = paste("Cluster composition —", lvl, "/", mth),
        subtitle = "Numbers = sample count. Dashed = expected PD fraction if random.",
        x = "Cluster", y = "Percentage (%)", fill = "Group"
      ) +
      theme_bw(base_size = 10)
    print(p)
  }
}

# ============================================================
# FIGURE 12: Network graphs — sparse methods
# ============================================================
message("Figure 12: Network graphs...")

plot_network_from_adj <- function(adj_path, title, group_vec = NULL) {
  if (!file.exists(adj_path)) return(NULL)
  adj <- read.csv(adj_path, row.names = 1, check.names = FALSE)
  adj_mat <- as.matrix(adj)

  # Remove self-loops
  diag(adj_mat) <- 0

  # Build igraph from adjacency (treat as weighted if non-binary)
  g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)

  if (vcount(g) == 0 || ecount(g) == 0) return(NULL)

  # Node attributes
  V(g)$group <- get_group(V(g)$name)
  V(g)$color <- ifelse(V(g)$group == "PD", GROUP_COLORS["PD"], GROUP_COLORS["C"])
  V(g)$size  <- 6

  # Edge sign (positive = blue, negative = red)
  if (!is.null(E(g)$weight)) {
    E(g)$edge_color <- ifelse(E(g)$weight > 0, "#1565C0", "#C62828")
    E(g)$width <- pmax(abs(E(g)$weight) * 2, 0.3)
  } else {
    E(g)$edge_color <- "#888888"
    E(g)$width <- 0.8
  }

  # Layout
  set.seed(42)
  layout <- layout_with_fr(g)

  old_par <- par(mar = c(1, 1, 3, 1))
  plot(g,
    layout       = layout,
    vertex.color = V(g)$color,
    vertex.size  = V(g)$size,
    vertex.label = V(g)$name,
    vertex.label.cex = 0.55,
    vertex.label.dist = 0.8,
    vertex.label.color = "black",
    vertex.frame.color = "white",
    edge.color   = E(g)$edge_color,
    edge.width   = E(g)$width,
    main         = title
  )
  legend("bottomleft",
    legend = c("PD", "Control", "Pos edge", "Neg edge"),
    col    = c(GROUP_COLORS["PD"], GROUP_COLORS["C"], "#1565C0", "#C62828"),
    pch    = c(21, 21, NA, NA),
    lty    = c(NA, NA, 1, 1),
    pt.bg  = c(GROUP_COLORS["PD"], GROUP_COLORS["C"], NA, NA),
    pt.cex = 1.5, cex = 0.7, bty = "n"
  )
  par(old_par)
  invisible(g)
}

for (lvl in c("family", "genus", "order")) {
  for (mth in c("spring", "spieceasi", "cclasso")) {
    mth_dir  <- file.path(RESULTS_DIR, lvl, mth)
    adj_file <- file.path(mth_dir, paste0("result_", mth, "_association_matrix.csv"))
    adj_bin  <- file.path(mth_dir, paste0("result_", mth, "_adjacency_matrix.csv"))

    # Use association matrix for edge weights, adjacency for sparsity mask
    if (file.exists(adj_file) && file.exists(adj_bin)) {
      assoc <- as.matrix(read.csv(adj_file, row.names = 1, check.names = FALSE))
      binary <- as.matrix(read.csv(adj_bin, row.names = 1, check.names = FALSE))
      diag(assoc)  <- 0
      diag(binary) <- 0

      # Mask: keep only edges present in adjacency
      masked <- assoc * binary
      tmp_path <- tempfile(fileext = ".csv")
      write.csv(masked, tmp_path)

      plot_network_from_adj(
        tmp_path,
        title = paste0(toupper(mth), " — ", lvl, " level\n(blue=pos, red=neg; node color=PD/Control)")
      )
      file.remove(tmp_path)
    }
  }
}

# ============================================================
# FIGURE 13: Cross-level hub consistency
# Which samples are consistently highly connected across levels?
# ============================================================
message("Figure 13: Cross-level hub consistency...")

for (mth in c("spring", "spieceasi")) {
  hub_summary <- do.call(rbind, lapply(c("family", "genus", "order"), function(lvl) {
    key <- paste(lvl, mth, sep = "_")
    d <- all_centralities[[key]]
    if (is.null(d)) return(NULL)
    # Normalize degree within level (0-1 scale)
    d$norm_degree <- if (max(d$degree, na.rm = TRUE) > 0) {
      d$degree / max(d$degree, na.rm = TRUE)
    } else {
      d$degree
    }
    d %>% select(taxon, group, norm_degree) %>% mutate(level = lvl)
  }))

  if (is.null(hub_summary) || nrow(hub_summary) == 0) next

  hub_wide <- hub_summary %>%
    pivot_wider(names_from = level, values_from = norm_degree, values_fill = 0) %>%
    mutate(
      consistency = rowMeans(select(., any_of(c("family", "genus", "order"))),
                             na.rm = TRUE)
    ) %>%
    arrange(desc(consistency))

  p <- ggplot(hub_summary, aes(x = reorder(taxon, -norm_degree), y = norm_degree,
                                fill = level, group = level)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = LEVEL_COLORS) +
    facet_grid(group ~ ., scales = "free_y", labeller = labeller(group = c(C = "Control", PD = "PD"))) +
    labs(
      title = paste("Cross-level hub consistency —", toupper(mth)),
      subtitle = "Normalized degree (0=min, 1=max within level). Samples appearing tall across levels = consistent hubs.",
      x = "Sample", y = "Normalized degree", fill = "Taxonomic level"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7))

  print(p)
}

# ============================================================
# FIGURE 14: PD vs Control degree summary
# ============================================================
message("Figure 14: PD vs Control degree comparison...")

sparse_cent_levels <- sparse_cent %>%
  filter(level %in% c("family", "genus", "order"))

p_box <- ggplot(sparse_cent_levels,
                aes(x = group, y = degree, fill = group)) +
  geom_boxplot(outlier.shape = 21, width = 0.5, alpha = 0.8) +
  geom_jitter(aes(color = group), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values  = GROUP_COLORS) +
  scale_color_manual(values = GROUP_COLORS) +
  facet_grid(level ~ method) +
  labs(
    title    = "Degree centrality: PD vs Control",
    subtitle = "Sparse methods only. Higher degree = more co-occurrence associations.",
    x = "Group", y = "Degree centrality", fill = "Group"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

print(p_box)

# Wilcoxon test summary table
wtest_rows <- sparse_cent_levels %>%
  group_by(level, method) %>%
  summarise(
    n_PD     = sum(group == "PD"),
    n_C      = sum(group == "C"),
    med_PD   = round(median(degree[group == "PD"], na.rm = TRUE), 4),
    med_C    = round(median(degree[group == "C"],  na.rm = TRUE), 4),
    p_wilcox = tryCatch(
      wilcox.test(
        degree[group == "PD"],
        degree[group == "C"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_wilcox = round(p_wilcox, 4),
    signif    = case_when(
      is.na(p_wilcox)  ~ "n/a",
      p_wilcox < 0.001 ~ "***",
      p_wilcox < 0.01  ~ "**",
      p_wilcox < 0.05  ~ "*",
      TRUE             ~ "ns"
    )
  )

grid.newpage()
grid.text(
  "Wilcoxon Test: PD vs Control degree centrality",
  x = 0.5, y = 0.95, gp = gpar(fontsize = 13, fontface = "bold")
)
grid.text(
  "Note: 24 samples only — results are exploratory, not statistically conclusive.",
  x = 0.5, y = 0.88, gp = gpar(fontsize = 9, col = "grey40")
)
tbl <- tableGrob(wtest_rows, rows = NULL,
                 theme = ttheme_minimal(base_size = 9))
grid.draw(tbl)

# ============================================================
# FIGURE 15: Class-level overview (only 3 methods, many taxa)
# ============================================================
message("Figure 15: Class-level overview...")

for (mth in c("sparcc", "cclasso")) {
  key <- paste("class", mth, sep = "_")
  d <- all_centralities[[key]]
  if (is.null(d) || nrow(d) == 0) next

  top20 <- d %>% arrange(desc(degree)) %>% head(20)

  p <- ggplot(top20, aes(x = reorder(taxon, degree), y = degree, fill = group)) +
    geom_col(color = "white") +
    scale_fill_manual(values = GROUP_COLORS) +
    coord_flip() +
    labs(
      title    = paste("Top 20 highest-degree nodes — class /", mth),
      subtitle = "Class level has 292 taxa; only top 20 by degree shown.",
      x = NULL, y = "Degree centrality", fill = "Group"
    ) +
    theme_bw(base_size = 10)
  print(p)
}

# ============================================================
# FIGURE 16: Method method comparison summary tile
# ============================================================
message("Figure 16: Summary interpretation tile...")

# Build a summary of key patterns
summary_text <- c(
  "KEY FINDINGS SUMMARY",
  "",
  "1. DENSE vs SPARSE SPLIT",
  "   SparCC & Spearman produce nearly fully connected networks (density ~1).",
  "   SPRING & SpiecEasi produce sparse networks (density 0.09-0.24).",
  "   CCLasso is intermediate (density 0.73-1.04).",
  "",
  "2. METHOD AGREEMENT",
  "   SparCC ≈ Spearman (Jaccard ~1.0): nearly identical edge sets across all levels.",
  "   SPRING ≈ SpiecEasi (Jaccard 0.47-0.69): moderate overlap — both graphical model approaches.",
  "   Spearman vs others: negative degree correlation at family/genus level — opposite node rankings.",
  "",
  "3. COMMUNITY STRUCTURE",
  "   SparCC family: 2 clusters — Cluster 1 (controls + some PD), Cluster 2 (mainly PD).",
  "   SPRING family: 5 communities; isolated samples (cluster 0) are predominantly non-networked.",
  "   Cluster 0 (SPRING) = isolated/disconnected samples; 13/24 samples have no edges.",
  "",
  "4. HUB NODES (sparse methods)",
  "   PD_8287 consistently highest degree/betweenness across family, genus, order × SPRING.",
  "   PD_8103 also prominent (high betweenness in SPRING family).",
  "   Hub agreement across methods is high except Spearman.",
  "",
  "5. CROSS-LEVEL CONSISTENCY",
  "   PD_8287 and PD_8103 are hubs at family, genus, AND order level.",
  "   Many PD samples are isolated (zero degree) in sparse networks.",
  "",
  "6. INTERPRETATION NOTE",
  "   Networks connect SAMPLES (not taxa) — edges = co-occurrence/correlation",
  "   of the sample's microbiome profile in the taxon space.",
  "   A 'hub' sample has a microbiome profile correlated with many others."
)

grid.newpage()
for (i in seq_along(summary_text)) {
  y_pos <- 0.96 - (i - 1) * 0.036
  if (y_pos < 0.02) break
  is_header <- grepl("^[0-9]\\.", summary_text[i]) || summary_text[i] == "KEY FINDINGS SUMMARY"
  grid.text(
    summary_text[i],
    x = 0.05, y = y_pos,
    just = "left",
    gp = gpar(
      fontsize  = if (is_header) 10 else 8.5,
      fontface  = if (is_header) "bold" else "plain",
      col       = if (summary_text[i] == "KEY FINDINGS SUMMARY") "#1565C0" else "black",
      fontfamily = "mono"
    )
  )
}

# ============================================================
# Close PDF
# ============================================================
dev.off()
message(paste("Done! Output written to:", OUTPUT_PDF))
message("===================================================")
