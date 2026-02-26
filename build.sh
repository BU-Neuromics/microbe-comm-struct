#!/bin/bash

# Build script for microbe-comm-struct Docker container
#
# NOTE: The bioc2u base image is amd64-only. On ARM64 hosts (e.g. DGX Spark),
# QEMU binfmt support must be registered first (one-time, persists until reboot):
#
#   docker run --privileged --rm tonistiigi/binfmt --install x86_64
#
# After that, this script builds an amd64 image that runs under QEMU emulation.

set -e  # Exit on error

echo "Building microbe-comm-struct Docker image (platform: linux/amd64)..."
docker build \
    --platform linux/amd64 \
    --progress=plain \
    -t microbe-comm-struct:latest \
    . 2>&1 | tee build.log

echo ""
echo "Build complete!"
echo "Image: microbe-comm-struct:latest (linux/amd64)"
echo ""
echo "To run an R script, use: ./run.sh path/to/your_script.R"
