#!/bin/bash
# Run an R script inside the microbe-comm-struct Docker container.
# The script's directory is mounted to /workspace inside the container.
#
# Usage:
#   ./run.sh <script.R> [args...]
#
# Examples:
#   ./run.sh 01_aitchison_pca.R counts.csv metadata.csv
#   ./run.sh 03_fastspar.R counts.csv metadata.csv

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Usage: $0 <script.R> [args...]"
    echo ""
    echo "Examples:"
    echo "  $0 01_aitchison_pca.R counts.csv metadata.csv"
    echo "  $0 03_fastspar.R counts.csv metadata.csv"
    exit 1
fi

R_SCRIPT="$1"
shift

if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: script not found: $R_SCRIPT"
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "$R_SCRIPT")" && pwd)
SCRIPT_NAME=$(basename "$R_SCRIPT")

echo "Running: $R_SCRIPT $*"
docker run --rm \
    --platform linux/amd64 \
    -v "$SCRIPT_DIR:/workspace" \
    microbe-comm-struct:latest \
    Rscript "/workspace/$SCRIPT_NAME" "$@"
