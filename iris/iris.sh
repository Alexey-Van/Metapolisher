#!/bin/bash
set -euo pipefail

if [ "$#" -lt 6 ]; then
  echo "Usage: $0 <input_dir> <output_dir> <genome.fa> <reads.bam> <vcf_in> <parliament_output>"
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
GENOME=$3
READS=$4
VCF_IN=$5
PARLIAMENT_OUTPUT=$6

IRIS_OUTPUT_DIR="IRIS_OUT"
VCF_OUT="sniffles_hifi_svs.filtered.iris.vcf"

# Download the IRIS image (if not already downloaded)
echo ">>> Pulling IRIS docker image..."
sudo docker pull quay.io/biocontainers/irissv:1.0.4--hdfd78af_2

echo ">>> Running IRIS refinement..."
sudo docker run -it \
-v "${INPUT_DIR}:${INPUT_DIR}" \
-v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
quay.io/biocontainers/irissv:1.0.4--hdfd78af_2 \
iris --hifi --keep_long_variants --keep_files \
genome_in="${GENOME}" \
reads_in="${READS}" \
vcf_in="${VCF_IN}" \
vcf_out="${OUTPUT_DIR}/${VCF_OUT}" \
out_dir="${OUTPUT_DIR}/${IRIS_OUTPUT_DIR}"

# Generating a list of SV files
echo ">>> Writing SV_filelist.txt..."
{
  echo "${PARLIAMENT_OUTPUT}"
  echo "${OUTPUT_DIR}/ONT_SVS_SNIFFLES_FILTERED"
  echo "${OUTPUT_DIR}/${VCF_OUT}"
} > "${OUTPUT_DIR}/SV_filelist.txt"

echo ">>> Done!"
