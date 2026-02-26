#!/usr/bin/env bash
# ============================================================
# run_pipeline.sh – launcher for the community structure pipeline
# ============================================================
# Runs all four analyses in parallel via Nextflow + Docker:
#   01  Aitchison PCA  (CLR-based ordination)
#   02  Robust PCA     (rCLR + softImpute + rsvd)
#   03  FastSpar       (co-occurrence network, bootstrap p-values)
#   04  Null-model     (Erdős–Rényi + configuration model validation)
#
# Prerequisites:
#   - nextflow  on PATH  (java >= 11 required)
#   - docker    on PATH  with microbe-comm-struct:latest built
#     (run ./build.sh if the image is missing)
#
#   On ARM64 hosts (e.g. DGX Spark), register QEMU before first use:
#     docker run --privileged --rm tonistiigi/binfmt --install x86_64
#
# Usage:
#   ./run_pipeline.sh --input <counts.csv> --metadata <meta.csv> [options]
#
# Required:
#   --input FILE       Counts CSV (rows = taxa, columns = samples)
#   --metadata FILE    Sample metadata CSV; must have a 'group' column
#
# Optional:
#   --taxonomy FILE    Taxonomy CSV for FastSpar heatmap annotation
#                      Columns: taxon_id, kingdom, phylum, class, order, family, genus
#   --output_dir DIR   Output directory [default: results]
#   --seed INT         Random seed [default: 42]
#
# Nextflow options:
#   -profile docker    Use Docker (default)
#   -profile local     Run without Docker (packages must be installed locally)
#   -profile slurm     Submit to SLURM via Singularity
#   -resume            Resume a previous run from cached task results
#
# Examples:
#   ./run_pipeline.sh --input counts.csv --metadata metadata.csv
#
#   ./run_pipeline.sh \
#       --input    counts.csv \
#       --metadata metadata.csv \
#       --taxonomy taxonomy.csv \
#       --output_dir results/family
#
#   ./run_pipeline.sh --input counts.csv --metadata metadata.csv -resume
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Validate nextflow is available ───────────────────────────
if ! command -v nextflow &>/dev/null; then
    echo "Error: nextflow not found on PATH."
    echo "Install via: conda install -c bioconda nextflow"
    echo "         or: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# ── Validate required arguments ──────────────────────────────
has_input=0; has_metadata=0
for arg in "$@"; do
    [[ "$arg" == "--input"    ]] && has_input=1
    [[ "$arg" == "--metadata" ]] && has_metadata=1
done

if (( ! has_input || ! has_metadata )); then
    grep "^#" "$0" | grep -v "^#!" | sed 's/^# \{0,1\}//' | head -50
    exit 1
fi

# ── Check Docker image ────────────────────────────────────────
if ! docker image inspect microbe-comm-struct:latest &>/dev/null; then
    echo "Docker image 'microbe-comm-struct:latest' not found — building now..."
    bash "$SCRIPT_DIR/build.sh"
fi

# ── Default to docker profile unless caller already passed -profile ──
profile_flag="-profile docker"
for arg in "$@"; do
    [[ "$arg" == "-profile" || "$arg" == "--profile" ]] && { profile_flag=""; break; }
done

echo ""
echo "============================================================="
echo "  Community Structure Analysis Pipeline"
echo "============================================================="
echo ""

nextflow run "$SCRIPT_DIR/main.nf" \
    $profile_flag \
    "$@"

echo ""
echo "Pipeline complete."
out_dir="$(echo "$@" | grep -oP '(?<=--output_dir\s)\S+' || true)"
out_dir="${out_dir:-results}"
echo "Results: ${out_dir}/"
