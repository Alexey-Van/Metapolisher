#!/bin/bash
set -euo pipefail

# Проверка аргументов
if [ "$#" -lt 3 ]; then
    echo "Использование: $0 genome.fa out.vcf vcf1.vcf vcf2.vcf ..."
    exit 1
fi

GENOME="$1"
OUTVCF="$2"
shift 2

# Создаём список VCF
VCF_LIST="vcfs.txt"
rm -f "$VCF_LIST"
for vcf in "$@"; do
    echo "$vcf" >> "$VCF_LIST"
done

# Запуск Jasmine
java -Xmx8g -jar jasmine.jar \
    "$VCF_LIST" \
    "$OUTVCF" \
    --genome "$GENOME" \
    --threads 8 \
    --output_genotypes
