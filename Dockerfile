FROM mambaorg/micromamba:1.5.8

USER root

# каналы
RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults

# ---------------------------
# mapping env
# ---------------------------
RUN micromamba create -y -n mapping \
    winnowmap \
    bwa-mem2 \
    samtools \
    && micromamba clean -a -y

# ---------------------------
# variant calling env
# ---------------------------
RUN micromamba create -y -n variant \
    sniffles=2.2 \
    bcftools \
    jasmine \
    && micromamba clean -a -y

# ---------------------------
# polishing env
# ---------------------------
RUN micromamba create -y -n polishing \
    merfin \
    repeatmasker \
    && micromamba clean -a -y

# ---------------------------
# python env
# ---------------------------
RUN micromamba create -y -n ml \
    python=3.11 \
    numpy \
    pandas \
    pysam \
    intervaltree \
    cyvcf2 \
    catboost \
    scikit-learn \
    && micromamba clean -a -y


# pipeline
COPY . /opt/pipeline/

RUN find /opt/pipeline -type f -name "*.sh" -exec cp {} /usr/local/bin/ \; && \
    chmod +x /usr/local/bin/*.sh

ENV PATH="/usr/local/bin:${PATH}"

USER mambauser