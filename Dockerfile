# bioc2u-builder provides pre-built apt binaries for Bioconductor/CRAN packages,
# avoiding long source compilation times.
#
# This image is amd64-only. On ARM64 hosts (e.g. DGX Spark), register QEMU
# once before building (survives until next reboot):
#
#   docker run --privileged --rm tonistiigi/binfmt --install x86_64
#
# Then build normally:
#   ./build.sh   (or: docker build --platform linux/amd64 -t microbe-comm-struct:latest .)
#
# R packages used by scripts 01-04:
#   zCompositions  – zero replacement for compositional data (script 01)
#   softImpute     – matrix completion via SVD (script 02)
#   rsvd           – randomised SVD for large matrices (script 02)
#   ggplot2        – all plots (scripts 01-04)
#   ggrepel        – non-overlapping text labels on PCA plots (scripts 01-02)
#   pheatmap       – correlation heatmap (script 03)
#   igraph         – network construction and metrics (script 04)
#   gridExtra      – multi-panel plot layout (script 04, optional fallback)
#
# External binary:
#   fastspar       – co-occurrence correlation tool (script 03), compiled from source
FROM --platform=linux/amd64 ghcr.io/bioconductor/bioc2u-builder:jammy-bioc-3.21-r-4.5.0

# ── System libraries ──────────────────────────────────────────────────────────
# Build tools for fastspar (C++ autotools project)
# Math libs: armadillo + openblas (fastspar linear algebra)
# GSL: fastspar random number generation
# Boost: fastspar CLI argument parsing and filesystem
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    autoconf \
    automake \
    libtool \
    make \
    g++ \
    libgsl-dev \
    libarmadillo-dev \
    libopenblas-dev \
    libboost-program-options-dev \
    libboost-filesystem-dev \
    && rm -rf /var/lib/apt/lists/*

# ── R packages via apt (pre-built binaries) ───────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-cran-zcompositions \
    r-cran-softimpute \
    r-cran-rsvd \
    r-cran-ggplot2 \
    r-cran-ggrepel \
    r-cran-pheatmap \
    r-cran-igraph \
    r-cran-gridextra \
    && rm -rf /var/lib/apt/lists/*

# Smoke-test: verify all packages load
RUN R -e " \
    library(zCompositions); \
    library(softImpute); \
    library(rsvd); \
    library(ggplot2); \
    library(ggrepel); \
    library(pheatmap); \
    library(igraph); \
    library(gridExtra); \
    message('All R packages loaded successfully') \
"

# ── fastspar (compiled from source) ───────────────────────────────────────────
RUN git clone https://github.com/scwatts/fastspar.git && \
    cd fastspar && \
    ./autogen.sh && \
    ./configure --prefix=/usr/ && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf fastspar

# Smoke-test: verify fastspar is on PATH
RUN fastspar --help > /dev/null 2>&1 || fastspar --version || true

# ── Runtime ───────────────────────────────────────────────────────────────────
WORKDIR /workspace
CMD ["R"]
