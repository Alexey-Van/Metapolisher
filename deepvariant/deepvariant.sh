#!/bin/bash

# Usage: ./deepvariant.sh <ref.fa> <input.bam> <out_dir> <model_type>
# Example: ./deepvariant.sh ref.fa sample.bam output_dir WGS 

set -e
set -o pipefail

ref=$1
bam=$2
out=$3
model=$4

if [[ -z $ref || -z $bam || -z $out || -z $model ]]; then
  echo "Usage: ./deepvariant.sh <ref.fa> <input.bam> <out_dir> <model_type>"
  echo "Model types: WGS, WES, PACBIO, ONT"
  exit 1
fi

mkdir -p $out

docker pull google/deepvariant:1.9.0

docker run --rm -v $PWD:/input -v $PWD/$out:/output \
  google/deepvariant:1.9.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=$model \
  --ref=/input/$ref \
  --reads=/input/$bam \
  --output_vcf=/output/output.vcf.gz \
  --output_gvcf=/output/output.g.vcf.gz \
  --num_shards=${SLURM_CPUS_PER_TASK:-8}
