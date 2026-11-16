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

RUN micromamba create -y -n sv \
    sniffles=2.8.1 \
    cutesv=2.1.4 \
    jasminesv=1.1.5 \
    samtools=1.25 \
    bcftools=1.25 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/sv/bin:$PATH

USER mambauser
