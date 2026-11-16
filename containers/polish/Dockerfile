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
    meryl=1.4.2 \
    yak=0.1 \
    k8=1.2 \
    merfin=1.1 \
    merqury=1.4.1 \
    samtools=1.25 \
    bcftools=1.25 \
    nextpolish2=0.2.2 \
    minimap2=2.31 \
    winnowmap=2.03 \
    racon=1.5.0 \
    medaka=2.2.2 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/polish/bin:$PATH

USER mambauser
