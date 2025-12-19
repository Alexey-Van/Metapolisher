#!/bin/bash

# Usage: ./parliament2.sh <input.bam> <input.bam.bai> <ref.fa> <ref.fa.fai> <out_dir>
# Example: ./parliament2.sh input.bam input.bam.bai chm13_chr1_mutated.fa chm13_chr1_mutated.fa.fai new

set -e
set -o pipefail

bam=$1
bai=$2
ref=$3
fai=$4
out=$5

if [[ -z $bam || -z $bai || -z $ref || -z $fai || -z $out ]]; then
  echo "Usage: ./parliament2.sh <input.bam> <input.bam.bai> <ref.fa> <ref.fa.fai> <out_dir>" 
  exit 1
fi

# Checking for file presence
for f in $bam $bai $ref $fai; do
  if [[ ! -s $f ]]; then
    echo "Error: file $f not found or empty"
    exit 1
  fi
done

mkdir -p $out

docker pull dnanexus/parliament2:0.1.11

docker run --rm \
  -v $PWD:/home/dnanexus/in \
  -v $PWD/$out:/home/dnanexus/out \
  dnanexus/parliament2:0.1.11 \
  --bam /home/dnanexus/in/$bam \
  --bai /home/dnanexus/in/$bai \
  -r /home/dnanexus/in/$ref \
  --fai /home/dnanexus/in/$fai \
  --prefix sample \
  --manta --lumpy --cnvnator \
  --delly_deletion --delly_insertion --delly_inversion --delly_duplication

echo ">>> Adding SOURCE=parliament to INFO fields..."

for vcf in $out/*.vcf; do
  if [[ -f "$vcf" ]]; then
    awk 'BEGIN{OFS="\t"} /^#/ {print; next} { $8=$8";SOURCE=parliament"; print }' "$vcf" > "${vcf}.tmp" && mv "${vcf}.tmp" "$vcf"
  fi
done