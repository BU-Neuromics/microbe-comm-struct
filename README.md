# microbe-comm-struct

Community structure analysis of microbial composition data derived from
ante-mortem cerebrospinal fluid (CSF) small RNA-seq in Parkinson's disease
patients and neurologically healthy controls.

The primary scientific question: do the observed microbial taxa show coherent
ecological community structure, or does the co-occurrence pattern resemble
random contamination?

---

## Methods

Four complementary analyses run in parallel via a Nextflow pipeline:

| Script | Method | Output |
|--------|--------|--------|
| `01_aitchison_pca.R` | Aitchison PCA (CLR-transformed PCA) | Scores + loadings PDFs |
| `02_rpca_deicode.R` | Robust PCA (rCLR + matrix completion) | Scores + loadings PDFs |
| `03_fastspar.R` | FastSpar co-occurrence network (500 iter, 1000 bootstrap) | Edge list CSV + heatmap PDF |
| `04_network_null_model.R` | Null-model validation vs Erdős–Rényi + configuration model | Comparison PDF |

> **Note:** An earlier approach using NetCoMi (SparCC, SPRING, SpiecEasi,
> CCLasso) was abandoned as all methods were computationally intractable on
> this dataset at every taxonomic level. See `archive/` for those scripts.

---

## Input data format

- **Counts CSV:** rows = samples, columns = taxa; row names = sample IDs,
  column names = taxon IDs. Non-integer (fractional) compositional counts.
- **Metadata CSV:** rows = samples; must have a `group` column
  (`PD` or `control`).
- **Taxonomy CSV (optional):** for FastSpar heatmap annotation.
  Columns: `taxon_id`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`.

---

## Quick start

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) with the `microbe-comm-struct:latest` image
- [Nextflow](https://www.nextflow.io/) (Java 11+)

**On ARM64 hosts (e.g. NVIDIA DGX Spark):** the Docker image is amd64-only.
Register QEMU once before building or running (persists until reboot):

```bash
docker run --privileged --rm tonistiigi/binfmt --install x86_64
```

### Build the Docker image

```bash
./build.sh
```

### Run the full pipeline

```bash
./run_pipeline.sh \
    --input    data/counts_family.csv \
    --metadata data/metadata.csv \
    --output_dir results/family
```

With optional taxonomy annotation for FastSpar:

```bash
./run_pipeline.sh \
    --input    data/counts_family.csv \
    --metadata data/metadata.csv \
    --taxonomy data/taxonomy.csv \
    --output_dir results/family
```

Resume after a failure (skips already-completed tasks):

```bash
./run_pipeline.sh --input counts.csv --metadata metadata.csv -resume
```

### Run a single script manually

```bash
./run.sh 01_aitchison_pca.R counts.csv metadata.csv
./run.sh 03_fastspar.R counts.csv metadata.csv
```

### Interactive shell in the container

```bash
./console.sh
```

---

## Pipeline outputs

```
results/
├── aitchison_pca/
│   ├── aitchison_scores.pdf      # PCA scores plot (colored by group, sized by library size)
│   └── aitchison_loadings.pdf    # Top 30 taxa by contribution to PC1+PC2
├── rpca_deicode/
│   ├── rpca_scores.pdf
│   └── rpca_loadings.pdf
├── fastspar/
│   ├── fastspar_edges.csv        # Filtered edge list (FDR < 0.05, |r| > 0.3)
│   └── fastspar_heatmap.pdf      # Correlation heatmap, top 50 taxa by degree
├── network_null/
│   └── null_model_comparison.pdf # Observed modularity/clustering vs null distributions
├── pipeline_report.html
├── pipeline_timeline.html
└── pipeline_trace.txt
```

---

## Nextflow profiles

| Profile | Description |
|---------|-------------|
| `docker` (default) | Runs all processes in `microbe-comm-struct:latest` |
| `local` | Skips Docker; packages must be installed locally |
| `slurm` | Submits to SLURM cluster via Singularity |

Configure resources in `nextflow.config`.

---

## Docker image

Base: `ghcr.io/bioconductor/bioc2u-builder:jammy-bioc-3.21-r-4.5.0` (amd64 only;
pre-built apt binaries for fast installs).

R packages: `zCompositions`, `softImpute`, `rsvd`, `ggplot2`, `ggrepel`,
`pheatmap`, `igraph`, `gridExtra`.

External binary: `fastspar` 1.0.0 (compiled from source).

---

## Repository layout

```
.
├── 01_aitchison_pca.R       # Aitchison PCA
├── 02_rpca_deicode.R        # Robust PCA / DEICODE
├── 03_fastspar.R            # FastSpar co-occurrence
├── 04_network_null_model.R  # Null-model validation
├── main.nf                  # Nextflow pipeline definition
├── nextflow.config          # Pipeline parameters and profiles
├── Dockerfile               # Container definition
├── build.sh                 # Build the Docker image
├── run.sh                   # Run a single script in the container
├── run_pipeline.sh          # Run the full Nextflow pipeline
├── console.sh               # Interactive shell in the container
├── environment.yml          # Conda env for Nextflow (openjdk + nextflow)
├── AGENTS.md                # Coding conventions for AI agents
└── archive/                 # Deprecated NetCoMi/SPRING scripts
```

---

## Running on a SLURM cluster

The Docker image is built automatically on every push to `main` via GitHub
Actions and published to GHCR:

```
ghcr.io/bu-neuromics/microbe-comm-struct:latest
```

Singularity (common on HPC) pulls from this registry automatically.

**If the package is private**, export credentials before running:

```bash
export SINGULARITY_DOCKER_USERNAME=<your_github_username>
export SINGULARITY_DOCKER_PASSWORD=<github_personal_access_token>
```

**Run on the cluster:**

```bash
./run_pipeline.sh \
    --input    data/counts_family.csv \
    --metadata data/metadata.csv \
    --output_dir results/family \
    -profile slurm
```

Adjust `queue`, `cpus`, `memory`, and `time` in the `slurm` profile of
`nextflow.config` to match your cluster's policies. For FastSpar with large
datasets, allocate at least 20 CPUs and 64 GB — the bootstrap loop
parallelizes across all available CPUs using `parallel::mclapply`.

> **Note:** The image is `linux/amd64` only. Standard x86_64 HPC nodes work
> natively. No QEMU needed on the cluster.

---

## CI/CD

A GitHub Actions workflow (`.github/workflows/docker-build.yml`) automatically:
- Builds the Docker image on every push to `main`
- Pushes `:latest` and `:sha-<short-sha>` tags to GHCR
- Creates `:v1.2.3` tags for version tags
- Uses GitHub Actions layer cache to speed up rebuilds

To pin a specific image version on the cluster instead of `:latest`, use the
`:sha-<short-sha>` tag (visible in the Actions run or `ghcr.io` package page).

---

## Citation

If using FastSpar, please cite:

- Watts SC et al. (2018) *FastSpar: rapid and scalable correlation estimation for compositional data.* Bioinformatics. doi:10.1093/bioinformatics/bty734
- Friedman J & Alm EJ (2012) *Inferring correlation networks from genomic survey data.* PLoS Comput Biol. doi:10.1371/journal.pcbi.1002687
