#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 5 ]; then
    echo "Usage: $0 <assembly.fa> <hifi.fq.gz> <short_R1.fq.gz> <short_R2.fq.gz> <threads>"
    exit 1
fi

ASM=$1
HIFI_BAM=$2
R1=$3
R2=$4
THREADS=$5

PREFIX=np2
LOG=nextpolish2.log

echo "[NP2] Assembly: $ASM" | tee $LOG
echo "[NP2] HiFi: $HIFI" | tee -a $LOG
echo "[NP2] Illumina: $R1 $R2" | tee -a $LOG
echo "[NP2] Threads: $THREADS" | tee -a $LOG

# 1) k‑меры из short‑reads → несколько yak
echo "[NP2] Building yak files from short reads" | tee -a $LOG

# Ты можешь менять k, b, split — здесь просто пример
yak count -b37 -t $THREADS -k21 -o ${PREFIX}.short.1.yak $R1 $R2 2>&1 | tee -a $LOG
yak count -b37 -t $THREADS -k31 -o ${PREFIX}.short.2.yak $R1 $R2 2>&1 | tee -a $LOG

# Собираем yak‑файлы в массив
YAKS=(
    ${PREFIX}.short.1.yak
    ${PREFIX}.short.2.yak
)

echo "[NP2] Yak files: ${YAKS[*]}" | tee -a $LOG

# 2) HiFi → assembly выравнивание
echo "[NP2] Mapping HiFi to assembly" | tee -a $LOG
minimap2 -t $THREADS -x map-hifi $ASM $HIFI \
  | samtools sort -@ $THREADS -o ${PREFIX}.hifi.bam 2>>$LOG
samtools index ${PREFIX}.hifi.bam

# 3) Запуск NextPolish2
echo "[NP2] Running nextPolish2" | tee -a $LOG
nextPolish2 \
  -t $THREADS \
  -o ${PREFIX}.polished.fa \
  $HIFI_BAM \
  $ASM \
  "${YAKS[@]}" 2>&1 | tee -a $LOG

echo "[NP2] Done. Output: ${PREFIX}.polished.fa" | tee -a $LOG
