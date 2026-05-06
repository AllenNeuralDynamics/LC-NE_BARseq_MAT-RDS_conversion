FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /
COPY r4-base.yml /etc/mamba-envs/r4-base.yml
COPY postInstall /
RUN /postInstall

# Python deps for the metadata-generation step (03_update_metadata.py).
# Installed against the base image's Python -- kept out of the r4-base conda
# env so the R and Python solver concerns stay independent.
RUN python3 -m pip install --no-cache-dir aind-data-schema==2.6.0
