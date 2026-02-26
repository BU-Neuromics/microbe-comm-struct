#!/usr/bin/env Rscript
# Usage: Rscript run_spring.R [ncores]
#   ncores - optional integer; number of CPU cores for SPRING StARS subsampling.
#            Defaults to max(1, detectCores() - 1).

# Install SPRING package if not already installed
if (!requireNamespace("SPRING", quietly = TRUE)) {
  cat("Installing SPRING package from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }
  remotes::install_github("GraceYoon/SPRING")
}

library(NetCoMi)
library(phyloseq)
library(SPRING)

# ── Core count ────────────────────────────────────────────────────────────────
# Accept an optional positional argument for the number of cores.
# When run via Nextflow, task.cpus is forwarded as args[1].
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && !is.na(suppressWarnings(as.integer(args[1])))) {
  n_cores <- max(1L, as.integer(args[1]))
} else {
  n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
message("Using ", n_cores, " core(s) for SPRING StARS subsampling.")

# Load data sets (same as run_netcomi.R)
# These are from the NetCoMi package
data("amgut1.filt") # ASV count matrix
data("amgut2.filt.phy") # phyloseq object

# Agglomerate to genus level
amgut_genus <- tax_glom(amgut2.filt.phy, taxrank = "Rank6")

# Rename taxonomic table and make Rank6 (genus) unique
amgut_genus_renamed <- renameTaxa(amgut_genus,
                                  pat = "<name>",
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Rank6")

head(tax_table(amgut_genus))

# Extract OTU table (count matrix) from phyloseq object
# SPRING expects n by p matrix where n = samples, p = OTUs
otu_matrix <- as(t(otu_table(amgut_genus_renamed)), "matrix")
# Ensure it's numeric
storage.mode(otu_matrix) <- "numeric"

# Print dimensions
cat("OTU matrix dimensions (samples x OTUs):", dim(otu_matrix), "\n")

# Filter to top 50 taxa by frequency (matching run_netcomi.R filtering)
# Calculate frequency (proportion of samples where OTU is present)
otu_freq <- colSums(otu_matrix > 0) / nrow(otu_matrix)
top50_idx <- order(otu_freq, decreasing = TRUE)[1:min(50, ncol(otu_matrix))]
otu_matrix_filtered <- otu_matrix[, top50_idx]

# Filter samples by total reads >= 1000 (matching run_netcomi.R filtering)
sample_totals <- rowSums(otu_matrix_filtered)
samples_keep <- sample_totals >= 1000
otu_matrix_filtered <- otu_matrix_filtered[samples_keep, ]

cat("Filtered matrix dimensions (samples x OTUs):", dim(otu_matrix_filtered), "\n")
cat("Number of samples kept:", sum(samples_keep), "\n")

# Run SPRING analysis
# Using quantitative = FALSE means compositional data (will apply mclr transformation)
# This is similar to the zero-handling approach in correlation-based network methods
cat("\nRunning SPRING analysis...\n")

fit_spring <- SPRING(data = otu_matrix_filtered,
                     quantitative = FALSE,  # Treat as compositional data
                     method = "mb",
                     lambda.min.ratio = 0.01,
                     nlambda = 20,
                     lambdaseq = exp(seq(log(0.6), log(0.6 * 0.01), length.out = 20)),
                     seed = 123456,  # Same seed as run_netcomi.R
                     ncores = n_cores,
                     thresh = 0.1,
                     subsample.ratio = 0.8,
                     rep.num = 20,
                     verbose = TRUE)

cat("\nSPRING analysis complete!\n")

# Display results summary
cat("\nOptimal lambda index:", fit_spring$output$stars$opt.index, "\n")
cat("Lambda sequence used:", length(fit_spring$lambdaseq), "values\n")

# Get the final network (precision matrix)
final_network <- fit_spring$fit$refit$stars

# Summary of the network
cat("\nNetwork summary:\n")
cat("Number of edges:", sum(final_network != 0) / 2, "\n")
cat("Network density:", sum(final_network != 0) / (ncol(final_network)^2 - ncol(final_network)), "\n")

# Save results
save(fit_spring, otu_matrix_filtered, file = "spring_results.RData")
cat("\nResults saved to spring_results.RData\n")
