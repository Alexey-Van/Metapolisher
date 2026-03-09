#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <input.bam> <input.fasta> <output_dir>" 
  exit 1
fi

INPUT_BAM="$1"
INPUT_FASTA="$2"
OUTPUT_DIR="$3"
OUTPUT_PREFIX="pepper"
THREADS=4

run_pepper_margin_deepvariant call_variant \
  -b "${INPUT_BAM}" \
  -f "${INPUT_FASTA}" \
  -o "${OUTPUT_DIR}" \
  -p "${OUTPUT_PREFIX}" \
  -t "${THREADS}" \
  --ont_r9_guppy5_sup
