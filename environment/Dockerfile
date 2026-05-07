FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /
COPY r4-base.yml /etc/mamba-envs/r4-base.yml

# Build the r4-base R conda environment.
# The `|| true` on the mamba update line lets the build continue if
# mamba's self-update fails -- the rest of the build doesn't depend on it.
RUN mamba update -n base -c conda-forge mamba libsolv -y || true \
 && rm -rf /opt/conda/envs/r4-base \
 && mamba env create -f /etc/mamba-envs/r4-base.yml \
 && echo "/opt/conda/envs/r4-base/lib" > /etc/ld.so.conf.d/r4-base.conf \
 && ldconfig

# Python deps for 02_update_metadata.py, installed against the base
# image's Python (separate from the r4-base conda env).
RUN python3 -m pip install --no-cache-dir aind-data-schema==2.7.1
