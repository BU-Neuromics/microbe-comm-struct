#!/usr/bin/env Rscript

# Script to run SPRING on any CSV file with taxa as columns and samples as rows

# Install and load required packages
if (!requireNamespace("SPRING", quietly = TRUE)) {
  message("Installing SPRING from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }
  remotes::install_github("GraceYoon/SPRING")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
  message("Installing optparse package...")
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}

library(SPRING)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to CSV file with taxa as columns and samples as rows [required]",
              metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="spring_output",
              help="Prefix for output files [default: %default]",
              metavar="PREFIX"),
  make_option(c("-c", "--ncores"), type="integer", default=NULL,
              help="Number of CPU cores to use [default: auto-detect and use N-1]",
              metavar="N"),
  make_option(c("-p", "--prevalence"), type="double", default=0.5,
              help="Minimum prevalence threshold (proportion of samples where taxon must be present) [default: %default]",
              metavar="PROP"),
  make_option(c("-t", "--threshold"), type="double", default=0.1,
              help="StARS threshold for stability selection [default: %default]",
              metavar="THRESH"),
  make_option(c("-n", "--nlambda"), type="integer", default=20,
              help="Number of lambda values for regularization path [default: %default]",
              metavar="N"),
  make_option(c("-r", "--rep-num"), type="integer", default=20,
              help="Number of StARS replications [default: %default]",
              metavar="N"),
  make_option(c("-s", "--subsample-ratio"), type="double", default=0.8,
              help="StARS subsample ratio [default: %default]",
              metavar="RATIO"),
  make_option("--seed", type="integer", default=123456,
              help="Random seed for reproducibility [default: %default]",
              metavar="SEED")
)

# Parse command line arguments
opt_parser <- OptionParser(
  option_list=option_list,
  description="\nRun SPRING (Sparse Partial correlation) network inference on taxonomic count data.",
  epilogue="Example:\n  Rscript run_spring_from_csv.R -i counts.csv -o results -c 8 -p 0.3"
)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("\nError: Input file (-i/--input) is required", call.=FALSE)
}

# Handle ncores argument
if (is.null(opt$ncores)) {
  # Auto-detect and use N-1 cores
  n_cores <- parallel::detectCores()
  n_cores_to_use <- max(1, n_cores - 1)
} else {
  n_cores_to_use <- opt$ncores
  if (n_cores_to_use < 1) {
    stop("Error: ncores must be a positive integer", call.=FALSE)
  }
}

# Validate other parameters
if (opt$prevalence < 0 || opt$prevalence > 1) {
  stop("Error: prevalence must be between 0 and 1", call.=FALSE)
}
if (opt$threshold < 0 || opt$threshold > 1) {
  stop("Error: threshold must be between 0 and 1", call.=FALSE)
}
if (opt$`subsample-ratio` <= 0 || opt$`subsample-ratio` > 1) {
  stop("Error: subsample-ratio must be between 0 and 1", call.=FALSE)
}

# Set variables from parsed arguments
input_csv <- opt$input
output_prefix <- opt$output
prevalence_threshold <- opt$prevalence
stars_threshold <- opt$threshold
nlambda <- opt$nlambda
rep_num <- opt$`rep-num`
subsample_ratio <- opt$`subsample-ratio`
seed <- opt$seed

# Check if input file exists
if (!file.exists(input_csv)) {
  stop(paste("Error: Input file", input_csv, "not found!"))
}

message("===================================================")
message("SPRING Network Inference from CSV")
message("===================================================")
message(paste("Input file:", input_csv))
message(paste("Output prefix:", output_prefix))
message(paste("Prevalence threshold:", prevalence_threshold))
message(paste("StARS threshold:", stars_threshold))
message(paste("Number of lambda values:", nlambda))
message(paste("StARS replications:", rep_num))
message(paste("StARS subsample ratio:", subsample_ratio))
message(paste("Random seed:", seed))
message("")

# Read the CSV file
message("Reading CSV file...")
data_raw <- read.csv(input_csv, row.names = 1, check.names = FALSE)

# Display data dimensions
message(paste("Data dimensions:", nrow(data_raw), "samples x", ncol(data_raw), "taxa"))
message("")

# Convert to numeric matrix (SPRING requirement)
data_matrix <- as.matrix(data_raw)
mode(data_matrix) <- "numeric"

# Check for and handle missing values
if (any(is.na(data_matrix))) {
  n_missing <- sum(is.na(data_matrix))
  message(paste("Warning:", n_missing, "missing values detected. Replacing with 0."))
  data_matrix[is.na(data_matrix)] <- 0
}

# Check for negative values
if (any(data_matrix < 0)) {
  stop("Error: Negative values detected in the data. SPRING requires non-negative counts.")
}

# Display basic statistics before filtering
message("Data summary (before filtering):")
message(paste("  Min value:", min(data_matrix)))
message(paste("  Max value:", max(data_matrix)))
message(paste("  Mean value:", round(mean(data_matrix), 2)))
message(paste("  Total zeros:", sum(data_matrix == 0), 
              paste0("(", round(100 * sum(data_matrix == 0) / length(data_matrix), 1), "%)")))
message("")

# Filter low-prevalence taxa
message("Filtering taxa...")
prevalence <- colSums(data_matrix > 0) / nrow(data_matrix)

keep_taxa <- prevalence >= prevalence_threshold

n_taxa_before <- ncol(data_matrix)
n_taxa_after <- sum(keep_taxa)
n_taxa_removed <- n_taxa_before - n_taxa_after

message(paste("  Prevalence threshold:", prevalence_threshold, "(taxa must be present in >=", 
              100 * prevalence_threshold, "% of samples)"))
message(paste("  Taxa before filtering:", n_taxa_before))
message(paste("  Taxa after filtering:", n_taxa_after))
message(paste("  Taxa removed:", n_taxa_removed))

if (n_taxa_removed > 0) {
  # Show some examples of removed taxa
  removed_taxa <- names(prevalence)[!keep_taxa]
  if (length(removed_taxa) <= 5) {
    message(paste("  Removed taxa:", paste(removed_taxa, collapse=", ")))
  } else {
    message(paste("  Example removed taxa:", paste(head(removed_taxa, 5), collapse=", "), "..."))
  }
}
message("")

# Apply the filter
if (n_taxa_after == 0) {
  stop("Error: No taxa remain after filtering! Try lowering the prevalence threshold.")
}

data_matrix <- data_matrix[, keep_taxa]

# Display statistics after filtering
message("Data summary (after filtering):")
message(paste("  Dimensions:", nrow(data_matrix), "samples x", ncol(data_matrix), "taxa"))
message(paste("  Total zeros:", sum(data_matrix == 0), 
              paste0("(", round(100 * sum(data_matrix == 0) / length(data_matrix), 1), "%)")))
message("")

# Run SPRING
message("Running SPRING...")
message("Parameters:")
message("  - quantitative: FALSE (compositional data)")
message("  - method: mb (Meinshausen-Bühlmann)")
message(paste("  - nlambda:", nlambda))
message("  - lambda.min.ratio: 0.01")
message(paste("  - StARS: subsample.ratio =", subsample_ratio, ", rep.num =", rep_num))
message(paste("  - threshold:", stars_threshold))
message(paste("  - ncores:", n_cores_to_use, "(", parallel::detectCores(), "total available)"))
message(paste("  - seed:", seed))
message("")

# Set up lambda sequence
lambda_max <- 0.6
lambda_min <- lambda_max * 0.01
lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

# Run SPRING with StARS selection
set.seed(seed)
fit_spring <- SPRING(
  data = data_matrix,
  quantitative = FALSE,      # Compositional data (applies mclr transformation)
  method = "mb",              # Meinshausen-Bühlmann method
  lambda.min.ratio = 0.01,
  nlambda = nlambda,
  lambdaseq = lambda_seq,
  seed = seed,
  ncores = n_cores_to_use,    # Use multiple cores for parallel computation
  thresh = stars_threshold,   # StARS threshold
  subsample.ratio = subsample_ratio,  # StARS subsample ratio
  rep.num = rep_num,          # Number of StARS replications
  verbose = TRUE
)

message("")
message("SPRING completed successfully!")
message("")

# Extract the optimal network
if (!is.null(fit_spring$output$stars$opt.index)) {
  opt_index <- fit_spring$output$stars$opt.index
  message(paste("Optimal lambda index:", opt_index))
  
  # Get the network adjacency matrix (from precision matrix)
  prec_matrix <- fit_spring$fit$refit$stars
  
  # Convert precision matrix to adjacency (non-zero off-diagonal elements)
  adj_matrix <- (prec_matrix != 0) * 1
  diag(adj_matrix) <- 0  # Remove self-loops
  
  # Network statistics
  n_edges <- sum(adj_matrix) / 2  # Divide by 2 because matrix is symmetric
  n_nodes <- ncol(adj_matrix)
  density <- n_edges / (n_nodes * (n_nodes - 1) / 2)
  
  message("")
  message("Network statistics:")
  message(paste("  Number of nodes (taxa):", n_nodes))
  message(paste("  Number of edges:", n_edges))
  message(paste("  Network density:", round(density, 4)))
  
  # Node degree distribution
  node_degrees <- colSums(adj_matrix)
  message(paste("  Average degree:", round(mean(node_degrees), 2)))
  message(paste("  Max degree:", max(node_degrees)))
  message(paste("  Isolated nodes:", sum(node_degrees == 0)))
} else {
  warning("No optimal network found by StARS selection")
  prec_matrix <- NULL
  adj_matrix <- NULL
}

# Save results
message("")
message("Saving results...")

# 1. Save the full SPRING fit object as RData
rdata_file <- paste0(output_prefix, "_spring_fit.RData")
save(fit_spring, data_matrix, file = rdata_file)
message(paste("  - Full results saved to:", rdata_file))

# 2. Save the precision matrix (network)
if (!is.null(prec_matrix)) {
  prec_file <- paste0(output_prefix, "_precision_matrix.csv")
  write.csv(prec_matrix, file = prec_file)
  message(paste("  - Precision matrix saved to:", prec_file))
  
  # 3. Save the adjacency matrix
  adj_file <- paste0(output_prefix, "_adjacency_matrix.csv")
  write.csv(adj_matrix, file = adj_file)
  message(paste("  - Adjacency matrix saved to:", adj_file))
  
  # 4. Save edge list (for network visualization tools)
  edges <- which(adj_matrix == 1, arr.ind = TRUE)
  edges <- edges[edges[,1] < edges[,2], ]  # Keep only upper triangle
  if (nrow(edges) > 0) {
    edge_list <- data.frame(
      from = colnames(data_matrix)[edges[,1]],
      to = colnames(data_matrix)[edges[,2]],
      weight = prec_matrix[edges]
    )
    edge_file <- paste0(output_prefix, "_edge_list.csv")
    write.csv(edge_list, file = edge_file, row.names = FALSE)
    message(paste("  - Edge list saved to:", edge_file))
  }
  
  # 5. Save node attributes
  node_attrs <- data.frame(
    taxon = colnames(data_matrix),
    degree = node_degrees,
    mean_abundance = colMeans(data_matrix),
    prevalence = colSums(data_matrix > 0) / nrow(data_matrix)
  )
  node_file <- paste0(output_prefix, "_node_attributes.csv")
  write.csv(node_attrs, file = node_file, row.names = FALSE)
  message(paste("  - Node attributes saved to:", node_file))
}

message("")
message("===================================================")
message("SPRING analysis completed successfully!")
message("===================================================")
