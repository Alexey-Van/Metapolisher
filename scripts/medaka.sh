#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <assembly.fa> <ont_reads.fq.gz> <threads>"
    exit 1
fi

ASM=$1
ONT_READS=$2
THREADS=$3

PREFIX=medaka
OUTDIR=${PREFIX}_output
LOG=medaka.log

echo "[MEDAKA] Assembly: $ASM" | tee "$LOG"
echo "[MEDAKA] ONT reads: $ONT_READS" | tee -a "$LOG"
echo "[MEDAKA] Threads: $THREADS" | tee -a "$LOG"

echo "[MEDAKA] Running medaka_consensus" | tee -a "$LOG"

medaka_consensus \
    -i "$ONT_READS" \
    -d "$ASM" \
    -o "$OUTDIR" \
    -t "$THREADS" \
    2>&1 | tee -a "$LOG"

echo "[MEDAKA] Done." | tee -a "$LOG"
echo "[MEDAKA] Output directory: $OUTDIR" | tee -a "$LOG"