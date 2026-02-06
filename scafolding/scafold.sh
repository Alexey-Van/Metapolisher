#!/bin/bash

# Использование:
#   ./run_yahs.sh contigs.fa hic_R1.fastq.gz hic_R2.fastq.gz prefix

# if [ $# -lt 4 ]; then
#     echo "Использование: $0 <contigs.fa> <hic_R1.fastq.gz> <hic_R2.fastq.gz> <prefix>"
#     exit 1
# fi

# CONTIGS=$1
# HIC_R1=$2
# HIC_R2=$3
# PREFIX=$4

# echo "FASTA контигов (референс для Hi-C): $CONTIGS"
# echo "Hi-C R1: $HIC_R1"
# echo "Hi-C R2: $HIC_R2"
# echo "Префикс: $PREFIX"

# # 1. Индексация FASTA
# echo "[1] samtools faidx..."
# samtools faidx "$CONTIGS"

# # 2. Индекс для minimap2
# echo "[2] minimap2 index..."
# minimap2 -d "${CONTIGS}.mmi" "$CONTIGS"

# # 3. Выравнивание Hi-C ридов на контиги и создание BAM
# echo "[3] Выравнивание Hi-C на контиги (minimap2)..."
# minimap2 -t 16 -ax sr "${CONTIGS}.mmi" "$HIC_R1" "$HIC_R2" | \
#     samtools view -b -o "${PREFIX}.hic.raw.bam" -

# # 4. Сортировка BAM по координатам (для удобства)
# echo "[4] Сортировка BAM по координатам..."
# samtools sort -o "${PREFIX}.hic.sorted.bam" "${PREFIX}.hic.raw.bam"
# samtools index "${PREFIX}.hic.sorted.bam"

# 5. Сортировка по имени для YaHS
echo "[5] Сортировка BAM по имени..."
samtools sort -n -o scafold_flye.hic.name_sorted.bam scafold_flye.hic.raw.bam

# # 6. Маркировка дупликатов
# echo "[6] Маркировка дупликатов (bammarkduplicates2)..."
# bammarkduplicates2 \
#     I="${PREFIX}.hic.name_sorted.bam" \
#     O="${PREFIX}.hic.markdup.bam" \
#     M="${PREFIX}.dup_metrics.txt"

samtools index scafold_flye.hic.name_sorted.bam

# 7. Запуск YaHS
echo "[7] Запуск YaHS..."
yahs -o scafold_flye data/T2T/original_fasta/flye/flye.contigs.fasta scafold_flye.hic.name_sorted.bam

echo "Готово."
echo "Основные выходы:"
# echo "  ${PREFIX}_scaffolds_final.fa"
# echo "  ${PREFIX}_scaffolds_final.agp"
