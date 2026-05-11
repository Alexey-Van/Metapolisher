FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && \
    apt-get install -y \
        procps \
        gawk \
        sed \
        grep \
        coreutils \
        findutils \
        util-linux \
        pigz \
        gzip \
        bzip2 \
        curl \
        wget \
        git \
        && apt-get clean
        
RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults

RUN micromamba create -y -n align \
    bwa-mem2 \
    samtools \
    minimap2 \
    winnowmap \
    racon \
    bcftools \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/align/bin:$PATH

USER mambauser