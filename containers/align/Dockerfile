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
    python=3.11 \
    bwa-mem3=0.17.0 \
    samtools=1.25 \
    bcftools=1.25 \
    pysam=0.24.1 \
    meryl=1.4.2 \
    minimap2=2.31 \
    winnowmap=2.03 \
    racon=1.5.0 \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/align/bin:$PATH

USER mambauser
