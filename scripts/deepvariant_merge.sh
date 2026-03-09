#!/bin/bash
set -euo pipefail

# Check arguments
if [ "$#" -lt 5 ]; then
  echo "Usage: $0 <reference.fa> <bam1> <bam2> <output_dir> <threads>" 
  exit 1
fi

REF=$1
BAM1=$2
BAM2=$3
OUTPUT_DIR=$4
THREADS=$5

# Name of the merged BAM
MERGED_BAM="${OUTPUT_DIR}/merged_hifi_illumina_hybrid.bam"
OUTPUT_VCF="${OUTPUT_DIR}/deepvariant_merged_hifi_illumina.vcf.gz"

# 1. BAM Association
echo ">>> Merging BAM files..."
samtools merge -@${THREADS} "${MERGED_BAM}" "${BAM1}" "${BAM2}"

# 2. BAM indexation
echo ">>> Indexing merged BAM..."
samtools index -@${THREADS} "${MERGED_BAM}"

# 3. Launch DeepVariant
docker pull google/deepvariant:1.9.0

echo ">>> Running DeepVariant..."
sudo docker run \
-v "$(dirname ${REF}):$(dirname ${REF})" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
google/deepvariant:1.9.0 \
/opt/deepvariant/bin/run_deepvariant \
--model_type HYBRID_PACBIO_ILLUMINA \
--ref "${REF}" \
--reads "${MERGED_BAM}" \
--output_vcf "${OUTPUT_VCF}" \
--num_shards "${THREADS}" \
--regions chr1 \
--intermediate_results_dir "${OUTPUT_DIR}/intermediate_results_dir"

echo ">>> Adding SOURCE=deepvariant to INFO fields..."

# DeepVariant produces a compressed VCF (.vcf.gz), so we need to process it with zcat/bgzip
TMP_VCF="${OUTPUT_DIR}/deepvariant_tmp.vcf"

zcat "${OUTPUT_VCF}" | \
awk 'BEGIN{OFS="\t"} /^#/ {print; next} { $8=$8";SOURCE=deepvariant"; print }' \
> "${TMP_VCF}"

bgzip -c "${TMP_VCF}" > "${OUTPUT_VCF}"
rm "${TMP_VCF}"
tabix -p vcf "${OUTPUT_VCF}"

echo ">>> Done! VCF with SOURCE=deepvariant saved to ${OUTPUT_VCF}"