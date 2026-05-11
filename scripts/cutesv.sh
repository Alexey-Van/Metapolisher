#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 <clr|ont|hifi> <bam> <ref.fasta> <out.vcf>"
    exit 1
fi

TYPE="$1"
BAM="$2"
REF="$3"
OUT="$4"

# Default parameters
INS_BIAS=""
INS_RATIO=""
DEL_BIAS=""
DEL_RATIO=""

case "$TYPE" in
    clr)
        INS_BIAS=100
        INS_RATIO=0.3
        DEL_BIAS=200
        DEL_RATIO=0.5
        ;;
    hifi)
        INS_BIAS=1000
        INS_RATIO=0.9
        DEL_BIAS=1000
        DEL_RATIO=0.5
        ;;
    ont)
        INS_BIAS=100
        INS_RATIO=0.3
        DEL_BIAS=100
        DEL_RATIO=0.3
        ;;
    *)
        echo "Error: type must be one of: clr, ont, hifi"
        exit 1
        ;;
esac

echo "[*] Running CuteSV with profile: $TYPE"
echo "    INS: bias=$INS_BIAS ratio=$INS_RATIO"
echo "    DEL: bias=$DEL_BIAS ratio=$DEL_RATIO"

mkdir -p cutesv_tmp
chmod +777 cutesv_tmp

docker run --rm \
    -v "$(pwd)":/data \
    aokad/cutesv:2.0.0 \
    cuteSV \
        /data/"$BAM" \
        /data/"$REF" \
        /data/"$OUT" \
        /data/cutesv_tmp \
        --max_cluster_bias_INS "$INS_BIAS" \
        --diff_ratio_merging_INS "$INS_RATIO" \
        --max_cluster_bias_DEL "$DEL_BIAS" \
        --diff_ratio_merging_DEL "$DEL_RATIO"

rm -rf cutesv_tmp
echo "[*] Done. Output: $OUT"
