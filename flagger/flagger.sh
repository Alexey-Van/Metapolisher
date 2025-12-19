#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <work_dir> <fasta_file> <bam_file>" 
  exit 1
fi

WORK_DIR="$1"
FASTA_FILE="$2"
BAM_FILE="$3"

cd "${WORK_DIR}"

# 1. Create FASTA index if missing
if [ ! -f "${FASTA_FILE}.fai" ]; then
  samtools faidx "${FASTA_FILE}"
else
  echo ">>> FAI index already exists: ${FASTA_FILE}.fai"
fi

# 2. Create BAM index if missing
if [ ! -f "${BAM_FILE}.bai" ]; then
  samtools index "${BAM_FILE}"
fi

# 3. Whole-genome BED
awk '{print $1"\t0\t"$2}' "${FASTA_FILE}.fai" > whole_genome.bed

# 4. JSON for annotations
echo "{" > annotations_path.json
echo "\"whole_genome\" : \"${PWD}/whole_genome.bed\"" >> annotations_path.json
echo "}" >> annotations_path.json

# 5. BAM → COV.GZ
docker run --rm -it \
  -v"${PWD}:${PWD}" mobinasri/flagger:v1.1.0--15d859f71ec26384837dee0add731b50aac158db \
  bam2cov \
    --bam "${PWD}/${BAM_FILE}" \
    --output "${PWD}/coverage_file.cov.gz" \
    --annotationJson "${PWD}/annotations_path.json" \
    --threads 16 \
    --baselineAnnotation whole_genome

# 6. Run hmm_flagger
mkdir -p hmm_flagger_outputs

             
docker run --rm -it \
  -v"${PWD}:${PWD}" mobinasri/flagger:v1.1.0 \
  hmm_flagger \
    --input "${PWD}/coverage_file.cov.gz" \
    --outputDir "${PWD}/hmm_flagger_outputs" \
    --iterations 50 \
    --threads 8

echo ">>> Results saved in ${PWD}/hmm_flagger_outputs"
