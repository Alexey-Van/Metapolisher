FROM mambaorg/micromamba:1.5.8

# Создаём окружение tools и ставим Sniffles2
RUN micromamba install -y -n base -c bioconda -c conda-forge sniffles=2.2 && \
    micromamba clean -a -y

# Устанавливаем остальные инструменты
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    sniffles=2.2 \
    jasmine \
    merfin \
    repeatmasker \
    python=3.10 \
    numpy \
    pandas \
    pysam \
    intervaltree \
    cyvcf2 \
    catboost \
    scikit-learn \
    && micromamba clean -a -y


USER root

COPY . /opt/pipeline/
RUN find /opt/pipeline -type f -name "*.sh" -exec cp {} /usr/local/bin/ \; && \
    chmod +x /usr/local/bin/*.sh

USER mambauser


ENV PATH="/usr/local/bin:${PATH}"
