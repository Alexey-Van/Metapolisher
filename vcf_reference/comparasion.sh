#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Usage example:
# ./comparasion.sh \
#   --truth data/T2T/original_fasta/truth.fixed.vcf \
#   --ml ml_filtered_9_trheshold_5_depth_10.vcf \
#   --fasta original.fasta \
#   --out results_dir
# -----------------------------

while [[ $# -gt 0 ]]; do
  case $1 in
    --truth) truth_vcf="$2"; shift 2 ;;
    --ml) ml_vcf="$2"; shift 2 ;;
    --fasta) fasta="$2"; shift 2 ;;
    --out) outdir="$2"; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done

mkdir -p "$outdir"

sort_and_index() {
    local infile="$1"
    local outfile="$2"

    echo "Sorting $infile → $outfile"
    bcftools sort "$infile" -Oz -o "$outfile"

    echo "Indexing $outfile"
    tabix -f -p vcf "$outfile"
}

# 1. Sort + index truth VCF
truth_gz="$outdir/truth.vcf.gz"
sort_and_index "$truth_vcf" "$truth_gz"

# 2. Sort + index ML VCF
ml_gz="$outdir/ml.vcf.gz"
sort_and_index "$ml_vcf" "$ml_gz"

# 3. bcftools isec
echo "Running bcftools isec..."
bcftools isec -p "$outdir/compare_out" "$truth_gz" "$ml_gz"

# # 4. Consensus
# echo "Building polished FASTA..."
# bcftools consensus -f "$fasta" "$ml_gz" > "$outdir/polished.fasta"

echo "Done. Results saved in: $outdir"
