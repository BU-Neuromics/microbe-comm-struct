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
# Nextflow is invoked via `conda run` using the environment defined in
# environment.yml (name: microbe-comm-struct).  The environment is created
# automatically on first use — no manual `conda activate` required.
#
# Prerequisites:
#   - conda (Miniconda or Anaconda) on PATH or in a standard location
#   - docker on PATH with microbe-comm-struct:latest built (or pulled from GHCR)
#
#   On ARM64 hosts (e.g. DGX Spark), register QEMU once before first use:
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
# Nextflow options (passed through directly):
#   -profile docker    Use Docker (default)
#   -profile local     Run without Docker (packages must be installed locally)
#   -profile sge       Submit to SGE cluster via Singularity
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
#   ./run_pipeline.sh --input counts.csv --metadata metadata.csv -profile sge
#   ./run_pipeline.sh --input counts.csv --metadata metadata.csv -resume
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONDA_ENV_NAME="microbe-comm-struct"
CONDA_ENV_FILE="$SCRIPT_DIR/environment.yml"

# ── Helpers ───────────────────────────────────────────────────
_info() { echo "[INFO]  $*"; }
_ok()   { echo "[OK]    $*"; }
_err()  { echo "[ERROR] $*" >&2; }

# ── 1. Locate conda ───────────────────────────────────────────
find_conda() {
    for candidate in \
        "$(command -v conda 2>/dev/null)" \
        "${CONDA_PREFIX:-}/bin/conda" \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/miniforge3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "/opt/conda/bin/conda"
    do
        if [[ -x "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

CONDA_CMD="$(find_conda)" || {
    _err "conda not found. Install Miniconda/Miniforge first, or add conda to PATH."
    exit 1
}
_ok "conda found: $CONDA_CMD"

# ── 2. Create the conda env if it doesn't exist ───────────────
if "$CONDA_CMD" env list | awk '{print $1}' | grep -qx "$CONDA_ENV_NAME"; then
    _ok "Conda env '$CONDA_ENV_NAME' already exists"
else
    _info "Conda env '$CONDA_ENV_NAME' not found — creating from $CONDA_ENV_FILE ..."
    "$CONDA_CMD" env create -f "$CONDA_ENV_FILE"
    _ok "Conda env '$CONDA_ENV_NAME' created"
fi

# ── 3. Parse arguments ────────────────────────────────────────
# Known pipeline args are captured into named variables; everything else is
# collected into NXF_ARGS and forwarded to nextflow as-is.
INPUT=""
METADATA=""
TAXONOMY=""
OUTPUT_DIR="results"
SEED=42
PROFILE=""
NXF_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)      INPUT="$2";      shift 2 ;;
        --metadata)   METADATA="$2";   shift 2 ;;
        --taxonomy)   TAXONOMY="$2";   shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --seed)       SEED="$2";       shift 2 ;;
        -profile|--profile)
                      PROFILE="$2";    shift 2 ;;
        *)            NXF_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$INPUT" || -z "$METADATA" ]]; then
    # Print only the contiguous comment block at the top of the file
    awk 'NR==1{next} /^[^#]/{exit} {sub(/^# ?/,""); print}' "$0"
    exit 1
fi

# ── 4. Check Docker image (docker profile only) ───────────────
if [[ "$PROFILE" != "sge" && "$PROFILE" != "local" ]]; then
    if command -v docker &>/dev/null; then
        if ! docker image inspect microbe-comm-struct:latest &>/dev/null; then
            _info "Docker image 'microbe-comm-struct:latest' not found — building..."
            bash "$SCRIPT_DIR/build.sh"
        else
            _ok "Docker image microbe-comm-struct:latest is present"
        fi
    fi
fi

# ── 5. Build profile flag ─────────────────────────────────────
if [[ -n "$PROFILE" ]]; then
    profile_flag="-profile $PROFILE"
else
    profile_flag="-profile docker"
fi

# ── 6. Verify Java works inside the conda env ─────────────────
# Resolve the env prefix so we can set NXF_JAVA_HOME explicitly.
# This sidesteps conflicts with system-wide JAVA_HOME on HPC clusters.
CONDA_ENV_PREFIX="$("$CONDA_CMD" run -n "$CONDA_ENV_NAME" bash -c 'echo $CONDA_PREFIX')"

# conda-forge openjdk can land in lib/jvm or expose java directly on PATH.
# Probe both locations and fall back to whatever JAVA_HOME the env exports.
if [[ -x "$CONDA_ENV_PREFIX/lib/jvm/bin/java" ]]; then
    CONDA_JAVA_HOME="$CONDA_ENV_PREFIX/lib/jvm"
elif [[ -x "$CONDA_ENV_PREFIX/bin/java" ]]; then
    CONDA_JAVA_HOME="$CONDA_ENV_PREFIX"
else
    # Last resort: ask the env what it thinks JAVA_HOME is
    CONDA_JAVA_HOME="$("$CONDA_CMD" run -n "$CONDA_ENV_NAME" bash -c 'echo "${JAVA_HOME:-}"')"
fi

if [[ -z "$CONDA_JAVA_HOME" || ! -x "$CONDA_JAVA_HOME/bin/java" ]]; then
    _err "Java not found inside conda env '$CONDA_ENV_NAME'."
    _err "Try removing and recreating the env:"
    _err "  conda env remove -n $CONDA_ENV_NAME && $0 --input ... --metadata ..."
    exit 1
fi

_ok "Java found: $("$CONDA_JAVA_HOME/bin/java" -version 2>&1 | head -1)"

# ── 7. Run pipeline via conda run ─────────────────────────────
echo ""
echo "============================================================="
echo "  Community Structure Analysis Pipeline"
echo "============================================================="
echo ""
_info "Launching Nextflow via conda env '$CONDA_ENV_NAME'..."
echo ""

# NXF_JAVA_HOME tells Nextflow exactly which Java to use, regardless of any
# system-wide JAVA_HOME set by the cluster environment or modules.
#
# Pipeline args are passed explicitly; NXF_ARGS carries any extra nextflow
# flags (e.g. -resume, -with-report) that were not consumed during parsing.
taxonomy_arg=""
[[ -n "$TAXONOMY" ]] && taxonomy_arg="--taxonomy $TAXONOMY"

"$CONDA_CMD" run -n "$CONDA_ENV_NAME" --no-capture-output \
    env NXF_JAVA_HOME="$CONDA_JAVA_HOME" \
    nextflow run "$SCRIPT_DIR/main.nf" \
    $profile_flag \
    --input      "$INPUT" \
    --metadata   "$METADATA" \
    --output_dir "$OUTPUT_DIR" \
    --seed       "$SEED" \
    $taxonomy_arg \
    "${NXF_ARGS[@]}"

echo ""
_ok "Pipeline complete. Results: ${OUTPUT_DIR}/"
