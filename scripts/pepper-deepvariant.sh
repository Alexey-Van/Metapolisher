#!/bin/bash
set -euo pipefail

INPUT_BAM="$1"
INPUT_FASTA="$2"
OUTPUT_DIR="$3"
OUTPUT_PREFIX="chm13_chr1" 
THREADS=4

CONTAINER_HOME="/home/pepper_run"
CONTAINER_BAM="$CONTAINER_HOME/ont_winnowmap.sorted.bam"
CONTAINER_FASTA="$CONTAINER_HOME/chm13_chr1_mutated.fa"
CONTAINER_OUT="$CONTAINER_HOME/output"

docker run --ipc=host \
  -v "$(dirname "$INPUT_BAM"):$CONTAINER_HOME" \
  -v "$(dirname "$INPUT_FASTA"):$CONTAINER_HOME" \
  -v "$OUTPUT_DIR:$OUTPUT_DIR" \
  kishwars/pepper_deepvariant:r0.8 \
  run_pepper_margin_deepvariant call_variant \
    -b "$CONTAINER_BAM" \
    -f "$CONTAINER_FASTA" \
    -o "$CONTAINER_OUT" \
    -p "$OUTPUT_PREFIX" \
    -t "$THREADS" \
    --ont_r9_guppy5_sup
