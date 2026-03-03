#!/usr/bin/env bash
set -e

### ============================================
### 1. СОЗДАНИЕ ОКРУЖЕНИЯ MICROMAMBA
### ============================================

# Создаём окружение
micromamba create -y -n sg_sandbox_env python=3.10
micromamba activate sg_sandbox_env

### ============================================
### 2. УСТАНОВКА ЗАВИСИМОСТЕЙ ЧЕРЕЗ MAMBA
### ============================================

micromamba install -y -c conda-forge -c bioconda \
    git \
    minimap2 \
    samtools \
    bedtools \
    seqtk \
    graphaligner \
    snakemake \
    biopython \
    parasail-python \
    pysam

### ============================================
### 3. КЛОНИРОВАНИЕ И СБОРКА sg_sandbox
### ============================================

git clone https://github.com/snurk/sg_sandbox
cd sg_sandbox

chmod +x build.sh
./build.sh

### ============================================
### 4. ПОДГОТОВКА ДАННЫХ
### ============================================

cd ..
mkdir -p data

if [ ! -f "data/hifi.fastq.gz" ]; then
    echo "Положите ваши HiFi риды в data/hifi.fastq.gz"
    exit 1
fi

# Тримминг HiFi ридов до 25 kbp
seqtk seq -L 25000 data/hifi.fastq.gz > data/hifi.trimmed.fastq

### ============================================
### 5. КОНФИГ sg_sandbox
### ============================================

cd sg_sandbox

cat > src/pipe/config.yaml <<EOF
S_MAX_RL: 20000
S_EXPECTED_COV: 30
EOF

### ============================================
### 6. ЗАПУСК ПОЛНОГО ПАЙПЛАЙНА STRING GRAPH
### ============================================

bash src/canu_launch/master_trim.sh \
    src/pipe/config.yaml \
    assembly_output \
    25000 \
    ../data/hifi.trimmed.fastq

### ============================================
### 7. РЕЗУЛЬТАТЫ
### ============================================

echo "Готово!"
echo "Основные файлы:"
echo "assembly_output/microasm/simplified.gfa"
echo "assembly_output/microasm/simplified.nodes.fasta"
