ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/c2-seurat:latest

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

COPY r4-base.yml /tmp/r4-base.yml

RUN mamba update -n base -c conda-forge mamba libsolv -y || true \
 && rm -rf /opt/conda/envs/* /etc/mamba-envs/* \
 && mamba env create -f /tmp/r4-base.yml \
 && rm /tmp/r4-base.yml \
 && ldconfig

RUN python3 -m pip install --no-cache-dir aind-data-schema==2.7.1
