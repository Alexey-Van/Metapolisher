#!/bin/bash

# Проверка аргумента
if [ -z "$1" ]; then
    echo "Использование: $0 <assembly.fasta>"
    exit 1
fi

ASSEMBLY="$1"
READS="$2"   
OUTDIR="lr_gapcloser_out"

mkdir -p "$OUTDIR"

LR_Gapcloser \
    -i "$ASSEMBLY" \
    -l "$READS" \
    -s ont \
    -o "$OUTDIR" \
    > "$OUTDIR/run.log" 2>&1
