#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <prefix> <ref.fa> <qry.fa> [threads]" >&2
    exit 2
fi

PREFIX="$1"
REF="$2"
QRY="$3"
THREADS="${4:-1}"

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer" >&2
    exit 2
fi

samtools faidx "$REF"

echo "[*] Aligning truth assembly to REF with minimap2..."
minimap2 -t "$THREADS" -ax asm5 --eqx "$REF" "$QRY" > "${PREFIX}.sam"
minimap2 -t "$THREADS" --cs -cx asm5 "$REF" "$QRY" > "${PREFIX}.paf"

echo "[*] Running SyRI..."
syri \
    -c "${PREFIX}.sam" \
    -r "$REF" \
    -q "$QRY" \
    --nc "$THREADS" \
    -F S \
    --prefix "$PREFIX"

SYRI_VCF="${PREFIX}syri.vcf"
if [ ! -s "$SYRI_VCF" ]; then
    echo "ERROR: SyRI did not create $SYRI_VCF" >&2
    exit 1
fi

PAFTOOLS="$(command -v paftools.js || true)"
if [ -z "$PAFTOOLS" ]; then
    echo "ERROR: paftools.js is not available on PATH" >&2
    exit 1
fi

echo "[*] Calling assembly variants with paftools..."
sort -k6,6 -k8,8n "${PREFIX}.paf" > "${PREFIX}.sorted.paf"
k8 "$PAFTOOLS" call -f "$REF" "${PREFIX}.sorted.paf" > "${PREFIX}.paftools.vcf"

echo "[*] Normalizing SyRI and paftools calls against REF..."
bcftools reheader --fai "${REF}.fai" \
    -o "${PREFIX}.syri.reheader.vcf" "$SYRI_VCF"

awk 'BEGIN { FS=OFS="\t" }
    /^#/ { print; next }
    $5 !~ /<[^>]+>/ { print }
' "${PREFIX}.syri.reheader.vcf" \
    | bcftools norm -f "$REF" -m -both -Ou \
    | bcftools sort -Oz -o "${PREFIX}.syri.sorted.vcf.gz"

bcftools norm -f "$REF" -m -both "${PREFIX}.paftools.vcf" -Ou \
    | bcftools sort -Oz -o "${PREFIX}.paftools.sorted.vcf.gz"

bcftools concat -a -d exact -Oz \
    -o "${PREFIX}.truth.concat.vcf.gz" \
    "${PREFIX}.syri.sorted.vcf.gz" \
    "${PREFIX}.paftools.sorted.vcf.gz"

bcftools sort -Oz \
    -o "${PREFIX}.truth.sorted.vcf.gz" \
    "${PREFIX}.truth.concat.vcf.gz"
bcftools norm -f "$REF" -d exact -Oz \
    -o truth.vcf.gz "${PREFIX}.truth.sorted.vcf.gz"
bcftools index -f truth.vcf.gz

echo "[*] Truth VCF created: truth.vcf.gz"
