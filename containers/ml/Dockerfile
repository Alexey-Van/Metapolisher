FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && \
    apt-get install -y \
        gcc \
        g++ \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        && rm -rf /var/lib/apt/lists/*

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

RUN micromamba create -y -n ml \
    python=3.11 \
    quast \
    numpy \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    scikit-learn \
    catboost \
    pysam \
    cyvcf2 \
    pyranges \
    intervaltree \
    pyarrow \
    joblib \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/ml/bin:$PATH

USER mambauser