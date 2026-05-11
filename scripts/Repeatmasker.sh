#!/bin/bash
set -euo pipefail

# ------------------------
# Check arguments
# ------------------------
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <repeatmasker_path> <genome.fa>"
  echo "  genome.fa         - input FASTA file"
  exit 1
fi

GENOME=$1


# ------------------------
# Run RepeatMasker
# ------------------------
RepeatMasker \
  -species human \
  -pa 8 \
  -dir zones/repeats \
  "$GENOME"

echo ">>> RepeatMasker finished. Results are in zones/repeats/"
