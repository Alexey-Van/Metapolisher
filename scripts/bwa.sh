#!/bin/bash

# Usage: ./bwa.sh <ref.fasta> <fastq_map> <out_prefix> [idx]
# fastq_map: prefix fastq 
 
set -o pipefail
set -e

ref=$1
fastq_map=$2
out=$3
wd=$PWD

if [[ -z $ref || -z $fastq_map || -z $out ]]; then
    echo "Usage: ./bwa.sh <ref.fasta> <fastq_map> <out_prefix> [idx]"
    exit 1
fi

# module load bwa
# module load samtools/1.21
conda install -y -c bioconda bwa=0.7.19
conda install -y -c bioconda samtools=1.22.1

cpu=${SLURM_CPUS_PER_TASK:-8}
idx=${SLURM_ARRAY_TASK_ID:-$4}
tmp=/lscratch/${SLURM_JOB_ID:-$$}

if [[ -z "$idx" ]]; then
    echo "No idx found (provide SLURM_ARRAY_TASK_ID or 4th argument)"
    exit 1
fi

out=${out}.${idx}

line=$(sed -n ${idx}p $fastq_map)
r1=$(echo $line | awk '{print $1}')
r2=$(echo $line | awk '{print $2}')

mkdir -p $tmp

# Main pipeline
bwa mem -t $cpu $ref $r1 $r2 > $tmp/${out}.sam

samtools fixmate -m -@$cpu $tmp/${out}.sam $tmp/${out}.fix.bam
rm $tmp/${out}.sam

samtools sort -@$cpu -O bam -o $tmp/${out}.bam -T $tmp/${out}.tmp $tmp/${out}.fix.bam
rm $tmp/${out}.fix.bam
samtools index -@$cpu $tmp/${out}.bam

samtools markdup -r -@$cpu $tmp/${out}.bam $tmp/${out}.dedup.bam
samtools index -@$cpu $tmp/${out}.dedup.bam

# Filtering 
samtools view -@$cpu -F0x100 -hb --write-index -o $tmp/${out}.dedup.pri.bam $tmp/${out}.dedup.bam

mv $tmp/${out}.dedup.pri.* ${wd}/
touch ${out}.done

echo "Done!"
