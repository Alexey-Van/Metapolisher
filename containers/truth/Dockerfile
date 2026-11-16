FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && \
    apt-get install -y \
        procps \
        gawk \
        coreutils \
        findutils \
        gzip \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults

RUN micromamba create -y -n truth \
    python=3.11 \
    syri=1.8.2 \
    k8=1.2 \
    minimap2=2.31 \
    samtools=1.25 \
    bcftools=1.25 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/truth/bin:$PATH

USER mambauser
