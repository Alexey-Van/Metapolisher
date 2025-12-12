#!/bin/bash
set -euo pipefail

# ------------------------
# Check arguments
# ------------------------
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <genome.fa> <meryl_db> <vcf1.vcf> [vcf2.vcf ... vcfN.vcf]"
  echo "  genome.fa   - reference FASTA"
  echo "  meryl_db    - k-mer database (created with meryl count)"
  echo "  vcf*.vcf    - one or more VCF files to annotate"
  exit 1
fi

GENOME=$1
MERYL_DB=$2
shift 2   # shift away genome and meryl_db, leaving only VCFs

# ------------------------
# Install Merfin (via bioconda)
# ------------------------
echo ">>> Installing Merfin via bioconda..."
conda install -y -c bioconda merfin

# ------------------------
# Process each VCF
# ------------------------
for VCF in "$@"; do
  if [ ! -f "$VCF" ]; then
    echo "File $VCF not found"
    exit 1
  fi

  echo ">>> Running Merfin on $VCF"

  # Run Merfin to evaluate variants
  merfin \
    -vcf "$VCF" \
    -ref "$GENOME" \
    -mers "$MERYL_DB" \
    -output "${VCF%.vcf}.merfin"

  # Suppose Merfin produces "${VCF%.vcf}.merfin.vcf" with MERFIN annotations.
  # Merge MERFIN INFO field back into the original VCF:
  bcftools annotate \
    -a "${VCF%.vcf}.merfin.vcf" \
    -c INFO/MERFIN \
    -o "${VCF%.vcf}.with_merfin.vcf" \
    -O v "$VCF"

  echo ">>> Annotated VCF saved to ${VCF%.vcf}.with_merfin.vcf"
done

echo ">>> All VCFs processed."
