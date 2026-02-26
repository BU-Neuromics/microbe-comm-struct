# bioc2u-builder provides pre-built apt binaries for Bioconductor/CRAN packages.
# This image is amd64-only. On ARM64 hosts (e.g. DGX Spark), register QEMU first:
#   docker run --privileged --rm tonistiigi/binfmt --install x86_64
# Then build with: docker build --platform linux/amd64 -t microbe-comm-struct:latest .
FROM --platform=linux/amd64 ghcr.io/bioconductor/bioc2u-builder:jammy-bioc-3.21-r-4.5.0

# Install system libraries required by R packages
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Install Bioconductor packages (phyloseq pulls in many Bioc deps; NetCoMi depends on limma, propr, etc.
RUN R -e "BiocManager::install(c( \
    'phyloseq', \
    'limma' \
    ), ask=FALSE, update=FALSE)"
RUN R -e "library(phyloseq); library(limma);"

RUN R -e "devtools::install_github('tpq/propr')"
RUN R -e "library(propr)"

## Install pulsar from CRAN archive (removed from CRAN 2026-02-15)
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pulsar/pulsar_0.3.11.tar.gz', repos=NULL, type='source')"
RUN R -e "library(pulsar)"

RUN R -e "install.packages(c('RcppEigen'), repos='https://cloud.r-project.org/')"
RUN R -e "library(RcppEigen)"

RUN R -e "devtools::install_github('zdk123/SpiecEasi')"
RUN R -e "library(SpiecEasi)"

RUN R -e "devtools::install_github('GraceYoon/SPRING')"
RUN R -e "library(SPRING)"

RUN R -e "devtools::install_github('stefpeschel/NetCoMi', \
            repos=c('https://cloud.r-project.org/', BiocManager::repositories()) \
         )"
RUN R -e "library(NetCoMi); message('All packages loaded successfully')"

RUN R -e "devtools::install_github('https://github.com/SarithaKodikara/LUPINE')"
RUN R -e "library(LUPINE)"

RUN R -e "install.packages(c('softImpute', 'zCompositions'), repos='https://cloud.r-project.org/')"

# for fastspar
RUN apt-get install -y libopenblas-dev

RUN git clone https://github.com/scwatts/fastspar.git; \
    cd fastspar; \
    ./autogen.sh; \
    ./configure --prefix=/usr/; \
    make; \
    make install

# Set default command to R
CMD ["R"]
