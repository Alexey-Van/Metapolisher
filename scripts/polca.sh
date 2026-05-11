#!/usr/bin/env bash
set -euo pipefail

# Использование:
#   ./polca.sh PREFIX ASM R1 R2
#
# Пример:
#   ./polca.sh HG002 HG002.asm.fa HG002.R1.fastq.gz HG002.R2.fastq.gz
#
# Ожидаемые файлы:
#   ASM — сборка (fasta)
#   R1  — парные риды R1
#   R2  — парные риды R2

if [ $# -lt 4 ]; then
  echo "Usage: $0 PREFIX ASM R1 R2"
  exit 1
fi

PREFIX="$1"
ASM="$2"
R1="$3"
R2="$4"

# Проверка входных файлов
for f in "$ASM" "$R1" "$R2"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: file not found or empty: $f" >&2
    exit 1
  fi
done

# Создаём выходную директорию
OUTDIR="pypolca_out"

# Запуск pypolca
# Опции взяты из официального CLI:
#   -a  сборка
#   -1  R1
#   -2  R2
#   -t  потоки
#   -o  выходная директория
#   -p  префикс
#   --careful  рекомендуемый режим (см. README)
pypolca run \
  -a "$ASM" \
  -1 "$R1" \
  -2 "$R2" \
  -t 8 \
  -o "$OUTDIR" \
  -p "$PREFIX" \
  --careful

echo "Done!"
echo "Corrected assembly: ${OUTDIR}/${PREFIX}_corrected.fasta"
echo "Report:             ${OUTDIR}/${PREFIX}.report"
