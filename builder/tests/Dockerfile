# Linux amd64 image for the fixture build (MMseqs2 13.45111 has no macOS arm64 build)
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8
# build context is builder/ (see Makefile)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN apt-get update && apt-get install -y --no-install-recommends make && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
WORKDIR /work
