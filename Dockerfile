# =============================================================================
# Dockerfile – H&E Spatial Heterogeneity Analysis
# Script: he_hypothesis_testing.R
# R: 4.1.1   |   Bioconductor: 3.14
# =============================================================================
#
# Build:
#   docker build -t he-analysis:latest .
#
# Run locally (bind-mount your data):
#   docker run --rm \
#     -v /path/to/h5files:/data/input:ro \
#     -v /path/to/results:/data/output \
#     he-analysis:latest \
#     --input_dir /data/input --output_dir /data/output
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_14

LABEL maintainer="sarah.andrews.park@gmail.com" \
      description="H&E tissue heterogeneity analysis pipeline" \
      r.version="4.1.1" \
      bioc.version="3.14"

# -----------------------------------------------------------------------------
# System libraries
#   libhdf5-dev   → rhdf5 / Rhdf5lib
#   libtiff-dev   → TIFF device used by ggsave(..., ".tiff")
#   libcurl4-openssl-dev / libssl-dev / libxml2-dev → common R package deps
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libtiff-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# remotes – needed to pin CRAN package versions
# -----------------------------------------------------------------------------
RUN Rscript -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"

# -----------------------------------------------------------------------------
# CRAN packages – pinned to Gerard's versions for reproducibility
# -----------------------------------------------------------------------------
RUN Rscript -e " \
  pkgs <- list( \
    c('dplyr',        '1.1.2'),   \
    c('tidyr',        '1.3.0'),   \
    c('ggplot2',      '3.5.1'),   \
    c('ggrepel',      '0.9.3'),   \
    c('scico',        '1.5.0'),   \
    c('tictoc',       '1.2.1'),   \
    c('Rfast',        '2.0.8'),   \
    c('nlme',         '3.1-168'), \
    c('compositions', '2.0-8')    \
  ); \
  for (p in pkgs) { \
    message('Installing ', p[1], ' ', p[2]); \
    remotes::install_version( \
      p[1], version = p[2], \
      repos   = 'https://cloud.r-project.org', \
      upgrade = 'never' \
    ) \
  }"

# -----------------------------------------------------------------------------
# Bioconductor packages – locked to BioC 3.14
# Install in dependency order: Rhdf5lib → rhdf5filters → rhdf5
# -----------------------------------------------------------------------------
RUN Rscript -e " \
  BiocManager::install( \
    c('Rhdf5lib', 'rhdf5filters', 'rhdf5'), \
    version = '3.14', \
    ask     = FALSE,  \
    update  = FALSE   \
  )"

# -----------------------------------------------------------------------------
# Copy analysis script
# -----------------------------------------------------------------------------
COPY R/he_hypothesis_testing.R /app/he_hypothesis_testing.R

WORKDIR /app

# -----------------------------------------------------------------------------
# Entry point
# Nextflow / AWS Batch will append CLI args, e.g.:
#   --input_dir /data/input --output_dir /data/output
# -----------------------------------------------------------------------------
ENTRYPOINT ["Rscript", "/app/he_hypothesis_testing.R"]
