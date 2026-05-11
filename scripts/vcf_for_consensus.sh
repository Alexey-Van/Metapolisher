#!/bin/bash
set -euo pipefail

# === Входные аргументы ===
VCF=$1
VCF_GZ="$1".gz
VCF_SORTED_GZ="$1".sorted.gz
VCF_SORTED_NORM_GZ="$1".norm.sorted.gz
DRAFT=$2
PREFIX=$3



gzip $VCF 
bcftools sort $VCF_GZ -Oz -o $VCF_SORTED_GZ
bcftools norm -m -both -f $DRAFT $VCF_SORTED_GZ -Oz -o $VCF_SORTED_NORM_GZ
bcftools index $VCF_SORTED_NORM_GZ
bcftools consensus -f $DRAFT $VCF_SORTED_NORM_GZ > $PREFIX.fasta