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

RUN micromamba create -y -n polish \
    meryl \
    merfin \
    merqury \
    samtools \
    bcftools \
    nextpolish2 \
    minimap2 \
    samtools \
    winnowmap \
    racon \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/polish/bin:$PATH

USER mambauser