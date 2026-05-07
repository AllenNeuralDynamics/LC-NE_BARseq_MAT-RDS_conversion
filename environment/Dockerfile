FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /
COPY r4-base.yml /etc/mamba-envs/r4-base.yml

# Build the r4-base R conda environment. mamba update is best-effort (it
# can fail under some build conditions and is non-critical). seurat-data
# is only published via GitHub, so it goes in via devtools after env
# creation, pinned to a specific commit for reproducibility.
RUN mamba update -n base -c conda-forge mamba libsolv -y || true \
 && rm -rf /opt/conda/envs/r4-base \
 && mamba env create -f /etc/mamba-envs/r4-base.yml \
 && echo "/opt/conda/envs/r4-base/lib" > /etc/ld.so.conf.d/r4-base.conf \
 && ldconfig \
 && Rscript -e '.libPaths("/opt/conda/envs/r4-base/lib/R/library/"); devtools::install_github("satijalab/seurat-data", ref = "3e51f44303069b64f5dc4d68e6a3d4a343f55c39")'

# Python deps for the metadata-generation step (02_update_metadata.py).
# Installed against the base image's Python -- kept out of the r4-base conda
# env so the R and Python environment builds stay independent.
RUN python3 -m pip install --no-cache-dir aind-data-schema==2.7.1
