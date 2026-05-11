#!/bin/bash
set -euo pipefail

# === Входные аргументы ===
DEEPVARIANT_VCF=$1
PEPPER_VCF=$2
DRAFT=$3
ILLUMINA_R1=$4
ILLUMINA_R2=$5
HIFI_READS=$6
OUT_PREFIX=$7

THREADS=${THREADS:-8}
K=${K:-21}

echo "[0] Starting pipeline..."
echo "Threads: $THREADS"
echo "k-mer size: $K"

# === 1. Фильтрация DeepVariant ===
echo "[1] Filtering DeepVariant..."
bcftools view -i 'FORMAT/VAF>=0.5 && FORMAT/GQ>=30' \
    "$DEEPVARIANT_VCF" -Oz -o ${OUT_PREFIX}.deepvariant.filtered.vcf.gz
bcftools index ${OUT_PREFIX}.deepvariant.filtered.vcf.gz

# === 2. Фильтрация PEPPER-DV ===
echo "[2] Filtering PEPPER-DV..."
bcftools view -i 'FORMAT/VAF>=0.5 && FORMAT/GQ>=25' \
    "$PEPPER_VCF" -Oz -o ${OUT_PREFIX}.pepper.filtered.vcf.gz
bcftools index ${OUT_PREFIX}.pepper.filtered.vcf.gz

# === 3. Объединение VCF ===
echo "[3] Merging filtered VCFs..."
bcftools norm -m-any -Oz -o ${OUT_PREFIX}.deepvariant.norm.vcf.gz ${OUT_PREFIX}.deepvariant.filtered.vcf.gz
bcftools norm -m-any -Oz -o ${OUT_PREFIX}.pepper.norm.vcf.gz ${OUT_PREFIX}.pepper.filtered.vcf.gz
bcftools index ${OUT_PREFIX}.deepvariant.norm.vcf.gz
bcftools index ${OUT_PREFIX}.pepper.norm.vcf.gz

bcftools merge \
    ${OUT_PREFIX}.deepvariant.norm.vcf.gz \
    ${OUT_PREFIX}.pepper.norm.vcf.gz \
    -Oz -o ${OUT_PREFIX}.merged.vcf.gz
bcftools index ${OUT_PREFIX}.merged.vcf.gz

# === 4. k‑меры сборки ===
echo "[4] Generating assembly k-mers..."
meryl count k=$K output ${OUT_PREFIX}.seqmers.meryl "$DRAFT"

# === 5. Illumina k‑меры (R1 + R2) ===
echo "[5] Generating Illumina k-mers..."
meryl count k=$K output ${OUT_PREFIX}.illumina_readmers.meryl \
    "$ILLUMINA_R1" "$ILLUMINA_R2"

# === 6. HiFi k‑меры ===
echo "[6] Generating HiFi k-mers..."
meryl count k=$K output ${OUT_PREFIX}.hifi_readmers.meryl "$HIFI_READS"

meryl union-sum \
    output ${OUT_PREFIX}.combined_readmers.meryl \
    ${OUT_PREFIX}.illumina_readmers.meryl \
    ${OUT_PREFIX}.hifi_readmers.meryl

# === 7. Merfin ===
echo "[7] Running Merfin..."
merfin -sequence $DRAFT -readmers ${OUT_PREFIX}.combined_readmers.meryl -seqmers ${OUT_PREFIX}.seqmers.meryl -peak 106 -output merfin_results -threads 8 -polish -vcf ${OUT_PREFIX}.merged.vcf.gz
