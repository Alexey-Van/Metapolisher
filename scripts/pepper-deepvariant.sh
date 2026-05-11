#!/bin/bash
set -euo pipefail

INPUT_BAM="$1"
INPUT_FASTA="$2"
OUTPUT_DIR="$3"
OUTPUT_PREFIX="$4"
THREADS="$5"

CONTAINER_HOME="/workspace"
CONTAINER_BAM="$CONTAINER_HOME/$(realpath --relative-to="$PWD" "$INPUT_BAM")"
CONTAINER_FASTA="$CONTAINER_HOME/$(realpath --relative-to="$PWD" "$INPUT_FASTA")"
CONTAINER_OUT="$CONTAINER_HOME/$OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

docker run --user root \
  -v "$PWD:$CONTAINER_HOME" \
  kishwars/pepper_deepvariant:r0.8 \
  run_pepper_margin_deepvariant call_variant \
    -b "$CONTAINER_BAM" \
    -f "$CONTAINER_FASTA" \
    -o "$CONTAINER_OUT" \
    -p "$OUTPUT_PREFIX" \
    -t "$THREADS" \
    --ont_r9_guppy5_sup
