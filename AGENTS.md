# Agent Guidelines for microbe-comm-struct

This document provides coding conventions and workflow guidelines for AI coding
agents working on the microbial community structure analysis codebase.

---

## Project Overview

Analysis of microbial composition data derived from **ante-mortem CSF small RNA
sequencing** in Parkinson's disease patients and neurologically healthy controls.
The scientific goal is community structure prediction: do the observed microbial
taxa show coherent ecological structure, or does the co-occurrence pattern resemble
random contamination?

### Brief history

The project started with the **NetCoMi R package**, which implements a broad family
of microbial co-occurrence network methods (SparCC, SPRING, CCLasso, etc.) and
provides a rich comparison infrastructure. However, all NetCoMi methods proved
**computationally intractable** at any taxonomic level for these data. Those scripts
now live in `archive/` for reference.

The current approach replaces NetCoMi with four **lighter-weight, standalone** R
scripts that are tractable on this dataset:

| Script | Method |
|--------|--------|
| `01_aitchison_pca.R` | Aitchison PCA (CLR-transformed PCA) |
| `02_rpca_deicode.R` | Robust PCA via rCLR + matrix completion (DEICODE-style) |
| `03_fastspar.R` | FastSpar co-occurrence network with bootstrap p-values |
| `04_network_null_model.R` | Null-model validation of the FastSpar network |

A **Nextflow pipeline** (`main.nf`) orchestrates all four scripts inside a Docker
container and produces HTML reports, timelines, and DAG diagrams.

---

## Repository Layout

```
.
├── 01_aitchison_pca.R         # Aitchison PCA (active)
├── 02_rpca_deicode.R          # Robust PCA / DEICODE (active)
├── 03_fastspar.R              # FastSpar co-occurrence (active)
├── 04_network_null_model.R    # Null-model validation (active)
├── main.nf                    # Nextflow pipeline
├── nextflow.config            # Pipeline parameters and profiles
├── run_pipeline.sh            # Convenience wrapper for nextflow
├── Dockerfile                 # Container definition
├── build.sh                   # Build the Docker image
├── run.sh                     # Run an R script inside the container
├── console.sh                 # Interactive R console in container
├── metadata.csv               # Sample metadata (group = PD / control)
├── environment.yml            # (Legacy) conda env spec — not actively used
├── NEW_METHODS.md             # Original specification for scripts 01-04
├── results/                   # Nextflow reports (HTML, trace, DAG)
├── work/                      # Nextflow task working directories
├── z2qm3/                     # OSF data mirror (counts CSVs, ANCOM-BC2 inputs)
└── archive/                   # DEPRECATED NetCoMi / SPRING scripts
    ├── run_netcomi.R
    ├── run_netcomi_method.R
    ├── run_spring.R
    ├── run_spring_from_csv.R
    └── run_lupine_single.R
```

---

## Build, Run & Test

### Build the Docker image
```bash
./build.sh
# Equivalent to: docker build --progress=plain -t microbe-comm-struct:latest .
```

### Run the full Nextflow pipeline (inside Docker)
```bash
./run_pipeline.sh \
  --input  z2qm3/osfstorage/ANCOM-BC2_count_input_files/filtered/filtered_ancom_count_input_genus.csv \
  --metadata metadata.csv \
  -profile docker
```

### Run a single R script manually
```bash
./run.sh 01_aitchison_pca.R counts.csv metadata.csv
./run.sh 03_fastspar.R counts.csv metadata.csv          # no taxonomy file
./run.sh 03_fastspar.R counts.csv metadata.csv taxa.csv # with taxonomy
```

### Interactive R console in the container
```bash
./console.sh
```

### Testing
No formal test suite. Validate by running individual scripts with sample data
from `z2qm3/` and inspecting PDF outputs and CSV edge lists.

---

## Script Descriptions

### `01_aitchison_pca.R`
- **Input:** counts CSV (rows = taxa, cols = samples), metadata CSV
- **Method:** Zero replacement via `zCompositions::cmultRepl` (CZM), then CLR
  transformation, then `prcomp`
- **Output:** `aitchison_scores.pdf`, `aitchison_loadings.pdf`
- **Stdout:** variance explained per PC (10 PCs), PC1/PC2 correlation with library size

### `02_rpca_deicode.R`
- **Input:** counts CSV, metadata CSV
- **Method:** rCLR (CLR using only non-zero entries per sample), then
  `softImpute` for matrix completion, then `rsvd::rsvd` for PCA
- **Output:** `rpca_scores.pdf`, `rpca_loadings.pdf`
- **Stdout:** chosen imputation rank, held-out reconstruction error

### `03_fastspar.R`
- **Input:** counts CSV, metadata CSV, (optional) taxonomy CSV
- **Method:** Calls the `fastspar` command-line binary; runs 500 iterations +
  1000 bootstraps for empirical p-values; FDR correction (BH)
- **Output:** `fastspar_edges.csv` (FDR < 0.05, |r| > 0.3), `fastspar_heatmap.pdf`
- **Stdout:** edge summary, positive/negative proportion, cross-kingdom breakdown

### `04_network_null_model.R`
- **Input:** edge list CSV from script 03, counts CSV
- **Method:** Builds an igraph network; computes modularity (Louvain), clustering
  coefficient, mean path length; compares to 1000 Erdős–Rényi and 1000
  configuration-model null graphs
- **Output:** `null_model_comparison.pdf`
- **Stdout:** z-scores and empirical p-values vs both null models, module sizes,
  positive/negative cohesion metrics

---

## Nextflow Pipeline Details

- **`main.nf`** defines four processes (AITCHISON_PCA, RPCA_DEICODE, FASTSPAR,
  NETWORK_NULL); NETWORK_NULL depends on the edge list from FASTSPAR.
- **`nextflow.config`** defines `docker`, `local`, and `slurm` profiles.
- The Docker container mounts `${projectDir}:/workspace` so scripts access data
  at `/workspace/...` inside the container.
- Optional taxonomy file: pass `--taxonomy path/to/taxonomy.csv`; if omitted the
  pipeline passes a dummy `NO_TAXONOMY` path to the process.
- Reports are written to `results/` (HTML report, timeline, trace, DAG).

---

## R Code Style

### Conventions
- **Indentation:** 2 spaces (no tabs)
- **Assignment:** `<-` (never `=`)
- **Pipe:** do NOT use `|>` or `%>%` — use explicit intermediate variables
- **Line length:** 80–100 chars (max 120)
- **Seed:** always `set.seed(42)` at the top of every script

### CLI arguments
Use positional `commandArgs(trailingOnly = TRUE)`; keep it simple.

### Messaging
- `message()` → informational to stderr
- `cat()` → formatted output to stdout
- Print section headers to make progress visible:
  ```r
  message("=== Aitchison PCA ===")
  ```

### File I/O
- Save results in multiple formats: `.RData` for R objects, `.csv` for
  tables/matrices, `.pdf` for plots.
- Use output filenames with clear prefixes matching the method name.

### Error handling
- `stop("Error: ...", call. = FALSE)` for fatal errors
- Check file existence before reading: `if (!file.exists(path)) stop(...)`
- Validate inputs early (dimensions, NA, negatives)

---

## Docker Environment

- **Base image:** `ghcr.io/bioconductor/bioc2u-builder:jammy-bioc-3.21-r-4.5.0`
  ⚠️ **ARM64 note:** This image is x86_64 only. On ARM64 hosts (e.g. DGX Spark)
  the build fails with "no match for platform in manifest". See the **Build
  Troubleshooting** section below.
- **Key R packages:**
  - Bioconductor: phyloseq, limma
  - GitHub: propr (tpq), SpiecEasi (zdk123), SPRING (GraceYoon), NetCoMi
    (stefpeschel), LUPINE (SarithaKodikara)
  - CRAN: pulsar (archive), RcppEigen, softImpute, zCompositions
  - System binary: `fastspar` (compiled from source)
- **Working dir inside container:** `/workspace`

### Build Troubleshooting — ARM64 (DGX Spark)

The `bioc2u-builder` base image is **amd64-only**. This machine is arm64.
We use QEMU emulation to build and run the amd64 image natively.

**One-time QEMU setup (survives until reboot):**
```bash
docker run --privileged --rm tonistiigi/binfmt --install x86_64
```
After this, `linux/amd64` appears in the Docker supported platforms list.

**Build and run (unchanged after QEMU setup):**
```bash
./build.sh          # now passes --platform linux/amd64 internally
./run.sh 01_aitchison_pca.R counts.csv metadata.csv
```

**⚠️ After a reboot**, re-run the QEMU registration command above before building or running the container.

**Performance note:** The Dockerfile installs several GitHub packages (propr, SpiecEasi, SPRING, NetCoMi, LUPINE) that must compile from source. Under QEMU emulation this is very slow. These packages are not needed by the active scripts (01–04). Removing them from the Dockerfile would drastically reduce build time — see "Slimming the Dockerfile" below.

**Slimming the Dockerfile (optional):**
Scripts 01–04 only need: `zCompositions`, `softImpute`, `rsvd`, `ggplot2`,
`pheatmap`, `igraph`, and the `fastspar` binary. All are apt-installable from
bioc2u except fastspar (C++ build). Removing NetCoMi/SPRING/SpiecEasi/LUPINE
would cut the build from potentially hours to ~10–20 minutes even under QEMU.

---

## Input Data

- **Counts CSVs:** `z2qm3/osfstorage/ANCOM-BC2_count_input_files/`
  - `filtered/` — prevalence-filtered versions (recommended for analysis)
  - `unfiltered/` — raw counts at each taxonomic level
  - Files available for levels: kingdom, phylum, class, order, family, genus, species, superkingdom
- **Metadata:** `metadata.csv` — rows = samples, must have a `group` column (PD / control)
- **Counts format:** rows = taxa, columns = samples; row names = taxon IDs; non-integer (fractional) compositional counts

---

## Common Pitfalls

1. Counts are **fractional** (non-integer) — do not assume integer counts
2. No prevalence/abundance filtering in scripts — filtering happens upstream
3. Do NOT pipe (`|>` or `%>%`); use explicit variables for debuggability
4. Always set `set.seed(42)` for reproducibility
5. Scripts output to the same directory as the input counts file by default
6. FastSpar writes intermediate files to a temp dir — clean up after runs
7. Archive scripts (`archive/`) are deprecated — do not modify or run them
8. The Docker image build fails on ARM64 — see Build Troubleshooting above

---

## Archive (Deprecated)

The `archive/` directory contains the original NetCoMi/SPRING/LUPINE scripts.
These are preserved for reference but are **not maintained**:
- `run_netcomi.R` — full NetCoMi pipeline with SparCC
- `run_netcomi_method.R` — parametric NetCoMi with multiple methods
- `run_spring.R` — SPRING with phyloseq input
- `run_spring_from_csv.R` — SPRING CLI tool from CSV
- `run_lupine_single.R` — LUPINE single-method run

Reason for deprecation: all methods were computationally intractable on the full
dataset even at the most coarse taxonomic levels.
