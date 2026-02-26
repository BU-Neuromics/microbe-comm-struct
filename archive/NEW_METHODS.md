**System / task prompt for coding agent:**

You are implementing a microbial co-occurrence and community structure analysis pipeline in R for a metatranscriptomic dataset. The dataset comes from small RNA sequencing of cerebrospinal fluid and uses a fractional count assignment scheme, so taxa counts are non-integer and compositionally structured. The primary scientific question is whether the observed microbial taxa show coherent community structure (evidence of real ecology) vs. random contamination.

**Context and assumptions:**
- All R package dependencies are already installed in the environment — do not install any packages
- Scripts should be standalone and runnable from the command line via `Rscript`
- Input is a counts matrix as a CSV file: rows are taxa, columns are samples, with row names as taxon IDs and column names as sample IDs
- A separate metadata CSV will also be provided: rows are samples (matching column names of counts matrix), with at minimum a column called `group` containing `PD` or `control`
- Counts may be non-integer (fractional) — do not assume integer counts
- The dataset has approximately 24 samples and potentially thousands to tens of thousands of taxa depending on taxonomic level
- Do not perform any prevalence or abundance filtering unless explicitly noted — the user will handle filtering decisions separately

---

**Script 1: `01_aitchison_pca.R`**

Implement Aitchison PCA (CLR-transformed PCA) with the following specifics:

- Accept two command-line arguments: path to counts CSV and path to metadata CSV
- Replace zeros using the `zCompositions` package with the `cmultRepl` function (method = `"CZM"`) before CLR transformation — do not use pseudocounts
- Apply CLR transformation manually (log of each value minus row mean of logs) after zero replacement, operating on the transposed matrix so that taxa are variables and samples are observations
- Run standard PCA using `prcomp` on the CLR-transformed sample × taxa matrix
- Produce two output plots saved as PDFs:
  - Scores plot (samples) colored by `group` from metadata, with sample IDs as labels, showing PC1 vs PC2 with variance explained in axis labels
  - Loadings plot showing the top 30 taxa by contribution to PC1 and PC2 combined (sum of squared loadings), as a labeled dot plot
- Also overlay library size as point size on the scores plot so the relationship between library size and ordination position is immediately visible
- Print to stdout: a table of variance explained per PC for the first 10 PCs, and Pearson correlations of PC1 and PC2 scores with library size

---

**Script 2: `02_rpca_deicode.R`**

Implement RPCA using the `deicode` approach via the `rlang` and matrix decomposition tools available in base R and the `rsvd` package:

- Accept two command-line arguments: path to counts CSV and path to metadata CSV
- Implement robust CLR (rCLR) transformation: for each sample, compute CLR using only the non-zero entries (ignore zeros in the log-mean calculation rather than imputing them)
- Apply matrix completion via `softImpute` package to fill in the zero positions after rCLR transformation, using the default rank determined by cross-validation (try ranks 1 through min(n,p)/2, select by minimizing held-out reconstruction error on a 20% random holdout of non-zero entries)
- Run PCA on the completed matrix using `rsvd::rsvd` for efficiency at large p
- Produce the same two plots as Script 1 (scores colored by group with library size as point size, top taxa loadings) saved as PDFs with `_rpca` suffix
- Print the chosen rank and held-out reconstruction error to stdout

---

**Script 3: `03_fastspar.R`**

Implement a FastSpar correlation analysis pipeline:

- Accept two command-line arguments: path to counts CSV and path to metadata CSV
- FastSpar is a command-line tool, not an R package — this script should prepare inputs, call FastSpar via `system()`, and parse outputs
- Write the counts matrix to a FastSpar-compatible TSV format (samples as rows, taxa as columns, with an `#OTU ID` header convention — check FastSpar documentation format)
- Call FastSpar with the following parameters: 500 iterations, 10 exclusion iterations, using a temporary directory for intermediates. Also call the FastSpar bootstrap pipeline to generate empirical p-values: generate 1000 bootstrap samples, compute correlations on each, then run the FastSpar p-value computation step
- Parse the output correlation matrix and p-value matrix back into R
- Apply FDR correction (BH method) to the p-value matrix
- Produce the following outputs:
  - A filtered edge list CSV containing only pairs with FDR < 0.05 and |correlation| > 0.3, with columns: taxon1, taxon2, correlation, pvalue, fdr
  - A correlation heatmap of the top 50 taxa by degree (number of significant edges) using `pheatmap`, clustered by complete linkage, saved as PDF
  - A summary printed to stdout: total significant edges, proportion positive vs negative, mean and SD of significant correlations, and a cross-kingdom breakdown if a taxonomy mapping file is provided as an optional third argument
- If a taxonomy CSV is provided as a third argument (columns: taxon_id, kingdom, phylum, class, order, family, genus), annotate the heatmap with kingdom-level color bars and print the cross-kingdom edge enrichment table

---

**Script 4: `04_network_null_model.R`**

Implement null model validation of the FastSpar network:

- Accept the filtered edge list CSV from Script 3 as the first argument, and the counts matrix CSV as the second argument
- Build an igraph network from the edge list, with edge weights as absolute correlation values and edge sign (positive/negative) as an edge attribute
- Compute the following observed network metrics: modularity (using `cluster_louvain`), clustering coefficient (global), mean path length (on largest connected component only), degree distribution, and proportion of edges that are negative
- Generate 1000 Erdős-Rényi null graphs with the same number of nodes and edges using `erdos.renyi.game`, and 1000 configuration model null graphs preserving the exact degree sequence using `sample_degseq` with the `vl` method
- For each null graph compute the same modularity and clustering coefficient
- Produce a plot with two panels: observed modularity vs null distribution (histogram with observed value marked), and observed clustering coefficient vs null distribution, saved as PDF
- Print to stdout: z-scores and empirical p-values for modularity and clustering coefficient against both null models, and the modularity, number of modules, and size of the three largest modules from the observed network
- Also compute and print network cohesion metrics: positive cohesion (mean of positive correlations weighted by their values) and negative cohesion (mean of absolute negative correlations), following the Herren & McMahon (2018) definition

---

**General coding standards for all scripts:**
- Use `argparse` or `commandArgs(trailingOnly=TRUE)` for argument parsing — keep it simple with positional arguments
- Write informative messages to stderr at each major step so progress is visible when running
- Save all outputs to the same directory as the input counts file unless a dedicated output directory is added as a future argument
- Use `ggplot2` for all plots except the pheatmap heatmap; apply a clean minimal theme
- Set a random seed of `42` at the top of every script
- Do not use the pipe operator `|>` or `%>%` — write explicit intermediate variables for clarity and debuggability
- Add a comment block at the top of each script describing inputs, outputs, and dependencies
