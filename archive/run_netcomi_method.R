#!/usr/bin/env Rscript

# Script to run NetCoMi network inference with different association methods
# Supports comparison of various association/correlation measures

# Install and load required packages
if (!requireNamespace("optparse", quietly = TRUE)) {
  message("Installing optparse package...")
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}

if (!requireNamespace("NetCoMi", quietly = TRUE)) {
  message("Installing NetCoMi package...")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org/")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  # Install dependencies first
  devtools::install_github("zdk123/SpiecEasi")
  devtools::install_github("GraceYoon/SPRING")
  devtools::install_github("stefpeschel/NetCoMi",
                           repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
}

library(NetCoMi)
library(phyloseq)
library(optparse)

# Patch for SpiecEasi::symBeta segfault with sparse beta matrices.
# SpiecEasi 1.1.2 calls t(beta) in C code which crashes on dgCMatrix objects
# returned by SPRING or SpiecEasi MB. We monkey-patch NetCoMi's internal
# .calcAssociation to wrap the beta in as.matrix() before symBeta is called.
#
# Two call sites are patched:
#   1. SPRING:     springres$output$est$beta[[opt.K]]
#   2. SpiecEasi:  SpiecEasi::getOptBeta(spiecres)
local({
  orig_fn   <- get(".calcAssociation", envir = asNamespace("NetCoMi"))
  orig_body <- body(orig_fn)

  # Walk the body AST and wrap known sparse-beta expressions with as.matrix().
  wrap_beta <- function(node) {
    if (is.call(node)) {
      node_str <- deparse(node)
      # SPRING beta (sparse dgCMatrix)
      if (identical(node_str, "springres$output$est$beta[[opt.K]]")) {
        return(call("as.matrix", node))
      }
      # SpiecEasi MB beta (sparse dgCMatrix from getOptBeta)
      if (identical(node_str, "SpiecEasi::getOptBeta(spiecres)")) {
        return(call("as.matrix", node))
      }
      node[] <- lapply(node, wrap_beta)
    }
    node
  }

  new_body <- wrap_beta(orig_body)
  body(orig_fn) <- new_body
  assignInNamespace(".calcAssociation", orig_fn, ns = "NetCoMi")
  message("  Applied SpiecEasi symBeta patch for sparse beta matrices.")
})

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input file (CSV or RData with phyloseq object) [required]",
              metavar="FILE"),
  make_option(c("-m", "--method"), type="character", default="sparcc",
              help="Association measure: pearson, spearman, bicor, sparcc, cclasso, ccrepe, spieceasi, spring, gcoda, propr [default: %default]",
              metavar="METHOD"),
  make_option(c("-o", "--output"), type="character", default="netcomi_output",
              help="Prefix for output files [default: %default]",
              metavar="PREFIX"),
  make_option(c("-t", "--taxrank"), type="character", default=NULL,
              help="Taxonomic rank to agglomerate (e.g., Rank6 for genus, Rank5 for family)",
              metavar="RANK"),
  make_option("--filttax", type="character", default="highestFreq",
              help="Taxa filtering method: highestFreq, highestVar, numbSamp, relAb [default: %default]",
              metavar="METHOD"),
  make_option("--filttaxpar", type="integer", default=50,
              help="Taxa filtering parameter (e.g., top N taxa) [default: %default]",
              metavar="N"),
  make_option("--filtsamp", type="character", default="totalReads",
              help="Sample filtering method: totalReads, numbTaxa, highestFreq [default: %default]",
              metavar="METHOD"),
  make_option("--filtsamppar", type="integer", default=1000,
              help="Sample filtering parameter (e.g., min total reads) [default: %default]",
              metavar="N"),
  make_option("--zeromethod", type="character", default="none",
              help="Zero handling: none, pseudo, pseudoZO, multRepl, alrEM, bayesMult [default: %default]",
              metavar="METHOD"),
  make_option("--normmethod", type="character", default="none",
              help="Normalization: none, TSS, CSS, COM, rarefy, VST, clr [default: %default]",
              metavar="METHOD"),
  make_option("--sparsmethod", type="character", default="none",
              help="Sparsification: none, t-test, bootstrap, threshold [default: %default]",
              metavar="METHOD"),
  make_option("--sparspar", type="double", default=0.3,
              help="Sparsification parameter (threshold or alpha) [default: %default]",
              metavar="VALUE"),
  make_option("--seed", type="integer", default=123456,
              help="Random seed for reproducibility [default: %default]",
              metavar="SEED"),
  make_option("--analyze", type="logical", default=TRUE,
              help="Run network analysis (centralities, clustering, etc.) [default: %default]",
              metavar="BOOL"),
  make_option("--clustmethod", type="character", default="cluster_fast_greedy",
              help="Clustering method: hierarchical, cluster_fast_greedy, cluster_louvain, etc. [default: %default]",
              metavar="METHOD"),
  make_option("--hubpar", type="character", default="eigenvector",
              help="Hub detection centrality: degree, betweenness, closeness, eigenvector [default: %default]",
              metavar="CENTR"),
  make_option("--verbose", type="integer", default=2,
              help="Verbosity level: 0 (silent), 1 (important), 2 (all) [default: %default]",
              metavar="LEVEL"),
  make_option(c("-c", "--ncores"), type="integer", default=NULL,
              help="Number of CPU cores for SPRING/SpiecEasi StARS subsampling. Defaults to max(1, detectCores() - 1).",
              metavar="N")
)

# Parse command line arguments
opt_parser <- OptionParser(
  option_list=option_list,
  description="\nRun NetCoMi network inference with different association methods.
  
Available association methods:
  - Correlation: pearson, spearman, bicor
  - Compositional: sparcc, cclasso, gcoda, propr
  - Network inference: ccrepe, spieceasi, spring",
  epilogue="Example:
  Rscript run_netcomi_method.R -i data.csv -m sparcc -o sparcc_results
  Rscript run_netcomi_method.R -i phyloseq.RData -m spring -t Rank6 -o spring_genus"
)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("\nError: Input file (-i/--input) is required", call.=FALSE)
}

# Check if input file exists
if (!file.exists(opt$input)) {
  stop(paste("Error: Input file", opt$input, "not found!"))
}

# Set variables from parsed arguments
input_file <- opt$input
method <- tolower(opt$method)
output_prefix <- opt$output
tax_rank <- opt$taxrank
filt_tax <- opt$filttax
filt_tax_par <- opt$filttaxpar
filt_samp <- opt$filtsamp
filt_samp_par <- opt$filtsamppar
zero_method <- opt$zeromethod
norm_method <- opt$normmethod
spars_method <- opt$sparsmethod
spars_par <- opt$sparspar
seed <- opt$seed
do_analyze <- opt$analyze
clust_method <- opt$clustmethod
hub_par <- opt$hubpar
verbose_level <- opt$verbose

# Resolve number of cores: CLI flag > auto-detect
if (!is.null(opt$ncores) && opt$ncores >= 1L) {
  n_cores <- as.integer(opt$ncores)
} else {
  n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
}

# Validate method
valid_methods <- c("pearson", "spearman", "bicor", "sparcc", "cclasso", 
                   "ccrepe", "spieceasi", "spring", "gcoda", "propr")
if (!(method %in% valid_methods)) {
  stop(paste("Error: Invalid method. Must be one of:", 
             paste(valid_methods, collapse=", ")))
}

message("===================================================")
message("NetCoMi Network Inference")
message("===================================================")
message(paste("Input file:", input_file))
message(paste("Association method:", method))
message(paste("Output prefix:", output_prefix))
message(paste("Random seed:", seed))
if (method %in% c("spring", "spieceasi")) {
  message(paste("CPU cores (StARS):", n_cores))
}
message("")

# Set random seed
set.seed(seed)

# Load data
message("Loading data...")
if (grepl("\\.csv$", input_file, ignore.case = TRUE)) {
  # Load CSV file
  count_data <- read.csv(input_file, row.names = 1, check.names = FALSE)
  
  # Create a simple phyloseq object (count matrix only).
  # The CSV is samples × taxa (rows = samples, cols = taxa), so pass
  # taxa_are_rows = FALSE.  NetCoMi's netConstruct computes associations
  # between the columns of the internal data matrix; with taxa in columns
  # this correctly produces a taxa × taxa association matrix.
  otu_mat <- otu_table(as.matrix(count_data), taxa_are_rows = FALSE)
  physeq <- phyloseq(otu_mat)
  
  message(paste("  Loaded CSV with", nrow(count_data), "samples and", 
                ncol(count_data), "taxa"))
  
} else if (grepl("\\.RData$|\\.rda$", input_file, ignore.case = TRUE)) {
  # Load RData file containing phyloseq object
  loaded_objects <- load(input_file)
  
  # Try to find phyloseq object
  physeq <- NULL
  for (obj_name in loaded_objects) {
    obj <- get(obj_name)
    if (inherits(obj, "phyloseq")) {
      physeq <- obj
      message(paste("  Loaded phyloseq object:", obj_name))
      break
    }
  }
  
  if (is.null(physeq)) {
    stop("Error: No phyloseq object found in RData file!")
  }
  
  message(paste("  Samples:", nsamples(physeq), "| Taxa:", ntaxa(physeq)))
  
} else {
  stop("Error: Input file must be CSV or RData format!")
}
message("")

# Taxonomic agglomeration if specified
if (!is.null(tax_rank)) {
  message(paste("Agglomerating to taxonomic rank:", tax_rank))
  
  if (is.null(tax_table(physeq))) {
    stop("Error: Taxonomic table not available for agglomeration!")
  }
  
  physeq <- tax_glom(physeq, taxrank = tax_rank)
  
  # Rename taxa to make them unique
  physeq <- renameTaxa(physeq,
                       pat = "<name>",
                       substPat = "<name>_<subst_name>(<subst_R>)",
                       numDupli = tax_rank)
  
  message(paste("  Taxa after agglomeration:", ntaxa(physeq)))
  message("")
}

# Prepare filtering parameters
filt_tax_par_list <- list()
if (filt_tax == "highestFreq") {
  filt_tax_par_list$highestFreq <- filt_tax_par
} else if (filt_tax == "highestVar") {
  filt_tax_par_list$highestVar <- filt_tax_par
} else if (filt_tax == "numbSamp") {
  filt_tax_par_list$numbSamp <- filt_tax_par
} else if (filt_tax == "relAb") {
  filt_tax_par_list$relAb <- filt_tax_par / 100  # Convert percentage to proportion
}

filt_samp_par_list <- list()
if (filt_samp == "totalReads") {
  filt_samp_par_list$totalReads <- filt_samp_par
} else if (filt_samp == "numbTaxa") {
  filt_samp_par_list$numbTaxa <- filt_samp_par
} else if (filt_samp == "highestFreq") {
  filt_samp_par_list$highestFreq <- filt_samp_par
}

# Prepare sparsification parameters
spars_par_list <- NULL
if (spars_method == "threshold") {
  spars_par_list <- list(threshold = spars_par)
} else if (spars_method == "t-test") {
  spars_par_list <- list(alpha = spars_par)
} else if (spars_method == "bootstrap") {
  spars_par_list <- list(alpha = spars_par)
}

# For MB-based methods (spring, spieceasi) the number of taxa fed to the
# graphical lasso must be strictly less than the number of samples, otherwise
# SpiecEasi::symBeta() receives a degenerate beta matrix and segfaults.
# Cap the taxa count to (nsamples - 1) when using these methods.
if (method %in% c("spring", "spieceasi") && filt_tax == "highestFreq") {
  n_samples <- nsamples(physeq)
  max_taxa <- n_samples - 1L
  if (filt_tax_par > max_taxa) {
    warning(paste0(
      "Method '", method, "' requires taxa < samples. ",
      "Reducing filttaxpar from ", filt_tax_par, " to ", max_taxa,
      " (nsamples - 1 = ", n_samples, " - 1)."
    ))
    filt_tax_par <- max_taxa
    filt_tax_par_list[[filt_tax]] <- filt_tax_par
  }
}

# Build measurePar for methods that need extra configuration.
# SPRING/SpiecEasi MB: StARS subsampling stability fails to select an
# optimal lambda on small datasets (n < 30) with default parameters,
# producing an empty opt.K and a NULL beta → symBeta() segfault.
# Tightening thresh and subsample.ratio resolves this for small n.
measure_par <- NULL
if (method %in% c("spring", "spieceasi")) {
  n_samples <- nsamples(physeq)
  if (n_samples < 30) {
    stars_thresh <- 0.05
    sub_ratio    <- 0.6
    message(paste0(
      "  Note: small dataset (n=", n_samples, ") – using ",
      "thresh=", stars_thresh, ", subsample.ratio=", sub_ratio,
      " for ", toupper(method), "/StARS."
    ))
  } else {
    stars_thresh <- 0.1
    sub_ratio    <- 0.8
  }
  if (method == "spieceasi") {
    # spiec.easi() accepts thresh/subsample.ratio inside pulsar.params, not
    # as top-level arguments.
    measure_par <- list(
      pulsar.params = list(
        thresh          = stars_thresh,
        subsample.ratio = sub_ratio,
        ncores          = n_cores
      )
    )
  } else {
    # SPRING passes these as top-level measurePar entries (handled by NetCoMi)
    measure_par <- list(
      thresh          = stars_thresh,
      subsample.ratio = sub_ratio,
      ncores          = n_cores
    )
  }
}

# Network construction
message("Constructing network...")
message("Parameters:")
message(paste("  - Measure:", method))
message(paste("  - Taxa filtering:", filt_tax, "=", filt_tax_par))
message(paste("  - Sample filtering:", filt_samp, "=", filt_samp_par))
message(paste("  - Zero method:", zero_method))
message(paste("  - Normalization:", norm_method))
message(paste("  - Sparsification:", spars_method, 
              if (!is.null(spars_par_list)) paste("(", names(spars_par_list)[1], "=", spars_par, ")", sep="") else ""))
message("")

net_construct <- netConstruct(
  data = physeq,
  taxRank = tax_rank,
  filtTax = filt_tax,
  filtTaxPar = filt_tax_par_list,
  filtSamp = filt_samp,
  filtSampPar = filt_samp_par_list,
  measure = method,
  measurePar = measure_par,
  zeroMethod = zero_method,
  normMethod = norm_method,
  sparsMethod = spars_method,
  thresh = if (spars_method == "threshold") spars_par else NULL,
  alpha = if (spars_method %in% c("t-test", "bootstrap")) spars_par else 0.05,
  adjust = if (spars_method %in% c("t-test", "bootstrap")) "adaptBH" else "none",
  dissFunc = "signed",
  verbose = verbose_level,
  seed = seed
)

message("")
message("Network construction completed!")
message("")

# Display network summary
message("Network summary:")
message(paste("  Number of nodes:", nrow(net_construct$countMat1)))
message(paste("  Number of edges:", sum(net_construct$adjaMat1 != 0) / 2))
message(paste("  Edge density:", 
              round(sum(net_construct$adjaMat1 != 0) / 
                    (nrow(net_construct$adjaMat1) * (nrow(net_construct$adjaMat1) - 1)), 4)))
if (!is.null(net_construct$assoMat1)) {
  pos_edges <- sum(net_construct$assoMat1[net_construct$adjaMat1 != 0] > 0, na.rm = TRUE)
  total_edges <- sum(net_construct$adjaMat1 != 0) / 2
  message(paste("  Positive edges:", pos_edges, "/", total_edges, 
                paste0("(", round(100 * pos_edges / total_edges, 1), "%)")))
}
message("")

# Network analysis
net_props <- NULL
if (do_analyze) {
  message("Analyzing network...")
  message("Parameters:")
  message(paste("  - Clustering method:", clust_method))
  message(paste("  - Hub detection centrality:", hub_par))
  message("")
  
  net_props <- netAnalyze(
    net = net_construct,
    centrLCC = TRUE,
    clustMethod = clust_method,
    hubPar = hub_par,
    hubQuant = 0.95,
    verbose = verbose_level
  )
  
  message("")
  message("Network analysis completed!")
  message("")
  
  # Display analysis summary
  message("Analysis summary:")
  
  # Component information
  comp_sizes <- net_props$compSize1
  message(paste("  Number of components:", ncol(comp_sizes)))
  message(paste("  Largest component size:", comp_sizes[1, 1]))
  
  # Global properties
  global_props <- net_props$globalPropsLCC
  if (!is.null(global_props$clustering) && is.numeric(global_props$clustering)) {
    message(paste("  LCC clustering coefficient:", round(global_props$clustering, 3)))
  }
  if (!is.null(global_props$modularity) && is.numeric(global_props$modularity)) {
    message(paste("  LCC modularity:", round(global_props$modularity, 3)))
  }
  if (!is.null(global_props$avPath) && is.numeric(global_props$avPath)) {
    message(paste("  LCC avg path length:", round(global_props$avPath, 3)))
  }
  
  # Hub information
  if (length(net_props$hubs$hubs1) > 0) {
    message(paste("  Number of hubs detected:", length(net_props$hubs$hubs1)))
    message(paste("  Hubs:", paste(net_props$hubs$hubs1, collapse=", ")))
  } else {
    message("  No hubs detected")
  }
  
  message("")
}

# Save results
message("Saving results...")

# 1. Save full results as RData
rdata_file <- paste0(output_prefix, "_", method, ".RData")
if (!is.null(net_props)) {
  save(net_construct, net_props, file = rdata_file)
} else {
  save(net_construct, file = rdata_file)
}
message(paste("  - Full results saved to:", rdata_file))

# 2. Save association matrix
if (!is.null(net_construct$assoMat1)) {
  asso_file <- paste0(output_prefix, "_", method, "_association_matrix.csv")
  write.csv(net_construct$assoMat1, file = asso_file)
  message(paste("  - Association matrix saved to:", asso_file))
}

# 3. Save adjacency matrix
adj_file <- paste0(output_prefix, "_", method, "_adjacency_matrix.csv")
write.csv(net_construct$adjaMat1, file = adj_file)
message(paste("  - Adjacency matrix saved to:", adj_file))

# 4. Save edge list
edges <- which(net_construct$adjaMat1 == 1, arr.ind = TRUE)
edges <- edges[edges[,1] < edges[,2], ]  # Keep only upper triangle
if (nrow(edges) > 0) {
  edge_list <- data.frame(
    from = rownames(net_construct$adjaMat1)[edges[,1]],
    to = colnames(net_construct$adjaMat1)[edges[,2]]
  )
  
  # Add association values if available
  if (!is.null(net_construct$assoMat1)) {
    edge_list$association <- net_construct$assoMat1[edges]
  }
  
  edge_file <- paste0(output_prefix, "_", method, "_edge_list.csv")
  write.csv(edge_list, file = edge_file, row.names = FALSE)
  message(paste("  - Edge list saved to:", edge_file))
}

# 5. Save network properties if analysis was performed
if (!is.null(net_props)) {
  # Save centralities
  # Build centrality data frame, handling possibly missing fields
  centr_list <- list(taxon = names(net_props$centralities$degree1))
  
  if (!is.null(net_props$centralities$degree1)) {
    centr_list$degree <- net_props$centralities$degree1
  }
  if (!is.null(net_props$centralities$between1)) {
    centr_list$betweenness <- net_props$centralities$between1
  }
  if (!is.null(net_props$centralities$close1)) {
    centr_list$closeness <- net_props$centralities$close1
  }
  if (!is.null(net_props$centralities$eigen1)) {
    centr_list$eigenvector <- net_props$centralities$eigen1
  }
  
  centr_df <- as.data.frame(centr_list)
  centr_file <- paste0(output_prefix, "_", method, "_centralities.csv")
  write.csv(centr_df, file = centr_file, row.names = FALSE)
  message(paste("  - Centralities saved to:", centr_file))
  
  # Save clusters
  if (!is.null(net_props$clustering$clust1)) {
    clust_df <- data.frame(
      taxon = names(net_props$clustering$clust1),
      cluster = net_props$clustering$clust1
    )
    clust_file <- paste0(output_prefix, "_", method, "_clusters.csv")
    write.csv(clust_df, file = clust_file, row.names = FALSE)
    message(paste("  - Clusters saved to:", clust_file))
  }
  
  # Save global properties
  if (!is.null(net_props$globalPropsLCC)) {
    global_df <- data.frame(
      property = names(net_props$globalPropsLCC),
      value = unlist(net_props$globalPropsLCC)
    )
    global_file <- paste0(output_prefix, "_", method, "_global_properties.csv")
    write.csv(global_df, file = global_file, row.names = FALSE)
    message(paste("  - Global properties saved to:", global_file))
  }
}

message("")
message("===================================================")
message("NetCoMi analysis completed successfully!")
message("===================================================")
