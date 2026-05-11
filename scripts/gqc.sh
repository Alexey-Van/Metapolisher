#!/usr/bin/env bash
set -euo pipefail

# Использование:
#   ./gqc_hap_compare.sh PREFIX query.fa ref.fa
#
# Пример:
#   ./gqc_hap_compare.sh CHM13 my_assembly.fa chm13v2.0.fa

if [ $# -lt 3 ]; then
    echo "Usage: $0 PREFIX QUERY_FASTA REF_FASTA"
    exit 1
fi

PREFIX="$1"
QUERY="$2"
REF="$3"

# Проверяем входные файлы
for f in "$QUERY" "$REF"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: file not found or empty: $f" >&2
        exit 1
    fi
done

# Создаём выходную директорию
mkdir -p "$PREFIX"

# Запуск GQC hapcompare (для гаплоидных сборок)
hapcompare \
    --qfasta "$QUERY" \
    --rfasta "$REF" \
    -p "$PREFIX" \
    -Q "query_${PREFIX}" \
    -R "ref_${PREFIX}"

echo "✔ GQC haploid comparison complete. Results saved in: $PREFIX/"
