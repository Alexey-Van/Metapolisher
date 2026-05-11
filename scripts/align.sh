#!/usr/bin/env bash
set -euo pipefail

# ==========================================
#   UNIVERSAL ALIGNER USING WINNOWMAP
#   Supports: hifi, clr, ont, illumina
#   Supports: .fastq, .fastq.gz, .fq, .fq.gz
#
#   Usage:
#     ./align.sh <type> <contigs.fasta> <reads...>
#
#   Examples:
#     ./align.sh hifi contigs.fasta pacbio.fastq.gz
#     ./align.sh ont contigs.fasta ont.fastq
#     ./align.sh illumina contigs.fasta R1.fq.gz R2.fq.gz
# ==========================================

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <hifi|clr|ont|illumina> <contigs.fasta> <reads...>"
    exit 1
fi

TYPE=$1
CONTIGS=$2
THREADS=8

mkdir -p ${TYPE}_bam

# ---------- Winnowmap requires a k-mer mask ----------
MASK="repetitive_k15.txt"

if [[ ! -f "$MASK" ]]; then
    echo "[INFO] Creating Winnowmap mask (repetitive k-mers)..."
    meryl count k=15 output meryl_db "$CONTIGS"
    meryl print greater-than 100 meryl_db > "$MASK"
fi

case "$TYPE" in

    hifi)
        READS=$3
        OUT="${TYPE}_bam/pacbio_hifi"
        echo "[INFO] Aligning PacBio HiFi with winnowmap..."
        winnowmap --MD -W "$MASK" -t $THREADS -ax map-pb "$CONTIGS" "$READS" > ${OUT}.sam
        ;;

    clr)
        READS=$3
        OUT="${TYPE}_bam/pacbio_clr"
        echo "[INFO] Aligning PacBio CLR with winnowmap..."
        winnowmap --MD -W "$MASK" -t $THREADS -ax map-pb "$CONTIGS" "$READS" > ${OUT}.sam
        ;;

    ont)
        READS=$3
        OUT="${TYPE}_bam/ont"
        echo "[INFO] Aligning ONT reads with winnowmap..."
        winnowmap --MD -W "$MASK" -t $THREADS -ax map-ont "$CONTIGS" "$READS" > ${OUT}.sam
        samtools view ${OUT}.sam -H > ${OUT}_header.sam
        samtools view ${OUT}.sam | awk 'length($6) < 64000' >> ${OUT}_header.sam
        OUT="${TYPE}_bam/ont_header"
        ;;

    illumina)
        R1=$3
        R2=$4
        OUT="${TYPE}_bam/illumina"
        echo "[INFO] Aligning Illumina PE reads with bwa-mem2..."
        bwa-mem2 index "$CONTIGS"
        bwa-mem2 mem -t $THREADS "$CONTIGS" "$R1" "$R2" > ${OUT}.sam
        ;;

    *)
        echo "Unknown type: $TYPE"
        exit 1
        ;;
esac

echo "[INFO] Converting SAM → sorted BAM..."
samtools sort ${OUT}.sam -@ $THREADS -o ${OUT}.sorted.bam
samtools index ${OUT}.sorted.bam

echo "[DONE] Output files:"
echo "  ${OUT}.sorted.bam"
echo "  ${OUT}.sorted.bam.bai"
