#!/bin/bash
set -euo pipefail

# ------------------------
# Check arguments
# ------------------------
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <repeatmasker_path> <genome.fa>"
  echo "  repeatmasker_path - directory where RepeatMasker is installed"
  echo "  genome.fa         - input FASTA file"
  exit 1
fi

REPEATMASKER_PATH=$1
GENOME=$2

# ------------------------
# Add RepeatMasker to PATH
# ------------------------
export PATH=$REPEATMASKER_PATH:$PATH

# ------------------------
# Run RepeatMasker
# ------------------------
RepeatMasker \
  -species human \
  -pa 8 \
  -dir zones/repeats \
  "$GENOME"

echo ">>> RepeatMasker finished. Results are in zones/repeats/"
