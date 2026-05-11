#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <work_dir> <fasta_file> <bam_file>" 
  exit 1
fi

FASTA_FILE="$1"
BAM_FILE="$2"

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
  -v"${PWD}:${PWD}" mobinasri/flagger:v1.1.0 \
  bam2cov \
    --bam "${PWD}/${BAM_FILE}" \
    --output "${PWD}/coverage_file.cov.gz" \
    --annotationJson "${PWD}/annotations_path.json" \
    --threads 4 \
    --baselineAnnotation whole_genome

# 6. Run hmm_flagger
mkdir -p hmm_flagger_outputs

             
docker run --rm   -v"${PWD}:/work" mobinasri/flagger:v1.1.0   hmm_flagger     --input /work/coverage_file.cov.gz     --outputDir /work/hmm_flagger_outputs --iterations 50     --threads 4

echo ">>> Results saved in ${PWD}/hmm_flagger_outputs"
