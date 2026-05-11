#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <prefix> <ref.fa> <qry.fa>"
    exit 1
fi

PREFIX="$1"

REF="$2"
QRY="$3"

OUTDIR="syri_${PREFIX}"
mkdir -p "${OUTDIR}"

echo "[*] Running minimap2..."
minimap2 -ax asm5 --eqx "${REF}" "${QRY}" > "${OUTDIR}/${PREFIX}.sam"

echo "[*] Running SyRI..."
syri \
  -c "${OUTDIR}/${PREFIX}.sam" \
  -r "${REF}" \
  -q "${QRY}" \
  --nc 8 -F S --prefix "${PREFIX}"

echo "[*] Done. Results in ${OUTDIR}"
