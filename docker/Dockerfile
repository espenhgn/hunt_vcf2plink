# hunt_vcf2plink: snakemake, plink2, triple-liftOver (with UCSC liftOver) and the workflow.
# Build from the repository root:
#   docker build --platform linux/amd64 -f docker/Dockerfile -t ghcr.io/espenhgn/hunt_vcf2plink:0.1.0 .
#
# UCSC liftOver and its chain files are free for non-commercial use only
# (https://genome.ucsc.edu/license/).

FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

USER root

# triple-liftOver has no releases; pin a commit. Its bundled liftOver binary is
# replaced by the pinned bioconda build.
ARG TRIPLE_LIFTOVER_COMMIT=f781ed6c5d016eac75b8be92938003207cdf5f0d
RUN python -c "import urllib.request; urllib.request.urlretrieve('https://github.com/GraceSheng/triple-liftOver/archive/${TRIPLE_LIFTOVER_COMMIT}.tar.gz', '/tmp/tlo.tar.gz')" && \
    mkdir -p /opt/triple-liftOver && \
    tar -xzf /tmp/tlo.tar.gz -C /opt/triple-liftOver --strip-components=1 && \
    rm /tmp/tlo.tar.gz && \
    ln -sf /opt/conda/bin/liftOver /opt/triple-liftOver/library/liftOver && \
    test -s /opt/triple-liftOver/library/chainfiles/hg19ToHg38.over.chain.gz && \
    chmod -R a+rX /opt/triple-liftOver

COPY workflow /opt/hunt_vcf2plink/workflow
COPY src/pca/pgen /opt/hunt_vcf2plink/src/pca/pgen
COPY config.yaml /opt/hunt_vcf2plink/config.yaml
RUN chmod -R a+rX /opt/hunt_vcf2plink

# apptainer exec does not run the micromamba entrypoint, so put the env on PATH
ENV PATH=/opt/conda/bin:$PATH

USER $MAMBA_USER
WORKDIR /tmp

LABEL org.opencontainers.image.source="https://github.com/espenhgn/hunt_vcf2plink" \
      org.opencontainers.image.description="HUNT VCF (GRCh37) to PLINK (GRCh38) and PCA workflow" \
      org.opencontainers.image.licenses="Apache-2.0"
