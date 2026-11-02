# hunt_vcf2plink: snakemake, plink2 and the workflow. UCSC liftOver and the chain
# file are downloaded by the workflow (rule setup_liftover), not part of the image.
# Build from the repository root:
#   docker build --platform linux/amd64 -f docker/Dockerfile -t ghcr.io/espenhgn/hunt_vcf2plink:0.1.0 .

FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

USER root

COPY workflow /opt/hunt_vcf2plink/workflow
COPY config.yaml /opt/hunt_vcf2plink/config.yaml
RUN chmod -R a+rX /opt/hunt_vcf2plink

# apptainer exec does not run the micromamba entrypoint, so put the env on PATH
ENV PATH=/opt/conda/bin:$PATH

USER $MAMBA_USER
WORKDIR /tmp

LABEL org.opencontainers.image.source="https://github.com/espenhgn/hunt_vcf2plink" \
      org.opencontainers.image.description="HUNT VCF (GRCh37) to PLINK (GRCh38) and PCA workflow" \
      org.opencontainers.image.licenses="Apache-2.0"
