#!/bin/bash
set -euo pipefail

# Check arguments
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <reference.fa> <input.bam> <output.vcf> <threads> <minsupport>"
  exit 1
fi

REF=$1
INPUT_BAM=$2
OUTPUT_VCF=$3
THREADS=$4
MINSUPPORT=$5

# Minimum length of SV
MINSVLEN=100

echo ">>> Installing Sniffles (if not already installed)..."
conda install -y -c bioconda sniffles=2.7.1

echo ">>> Running Sniffles..."
sniffles \
  --input "${INPUT_BAM}" \
  --vcf "${OUTPUT_VCF}" \
  --reference "${REF}" \
  --threads "${THREADS}" \
  --minsupport "${MINSUPPORT}" \
  --minsvlen "${MINSVLEN}"

echo ">>> Adding SOURCE=sniffles to INFO field..."
awk 'BEGIN{OFS="\t"} /^#/ {print; next} { $8=$8";SOURCE=sniffles"; print }' "${OUTPUT_VCF}" > "${OUTPUT_VCF}.tmp" && mv "${OUTPUT_VCF}.tmp" "${OUTPUT_VCF}"

echo ">>> Done!"
