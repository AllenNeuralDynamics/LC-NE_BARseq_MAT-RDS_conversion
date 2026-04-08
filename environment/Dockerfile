# hash:sha256:046cb60d5a2ac29373a12426c645a03516f6b237e844f5965082ee071ac6f252
FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /

COPY postInstall /
RUN /postInstall
