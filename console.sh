#!/bin/bash
# Open an interactive bash shell inside the microbe-comm-struct Docker container.
# The current directory is mounted to /workspace.
#
# Usage:
#   ./console.sh

set -euo pipefail

echo "Opening shell in microbe-comm-struct container..."
echo "Current directory mounted at /workspace"
echo ""
docker run --rm -it \
    --platform linux/amd64 \
    -v "$(pwd):/workspace" \
    microbe-comm-struct:latest \
    /bin/bash
