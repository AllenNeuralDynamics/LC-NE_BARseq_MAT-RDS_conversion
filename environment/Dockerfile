FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /
COPY r4-base.yml /etc/mamba-envs/r4-base.yml
COPY postInstall /
RUN /postInstall