#!/bin/bash

# Usage: ./winnowmap.sh <ref.fasta> <reads.fastq> <out_prefix> <map_mode> <extra_opts>
# Example: ./winnowmap.sh ref.fa reads.fq out_prefix map-pb ""

set -o pipefail
set -e

ref=$1
reads=$2
out=$3
map=$4        # example: map-pb or map-ont
opt=$5        # additional options (can be left blank)

cpus=${SLURM_CPUS_PER_TASK:-8}
wd=$PWD

if [[ -z $ref || -z $reads || -z $out || -z $map ]]; then
    echo "Usage: ./winnowmap.sh <ref.fasta> <reads.fastq> <out_prefix> <map_mode> <extra_opts>"
    exit 1
fi

# module load meryl
# module load winnowmap
# module load samtools/1.21

conda install -y -c bioconda meryl=1.4.1
conda install -y -c bioconda winnowmap=2.03
conda install -y -c bioconda samtools=1.22.1


# 1. Collecting repeating 15-mers
meryl count k=15 $ref output merylDB
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

# 2. Read alignment
winnowmap --MD -W repetitive_k15.txt -ax $map $opt -t$cpus $ref $reads > ${out}.sam

# 3. Sorting and converting to BAM
samtools sort -@$cpus -m2G -T ${out}.tmp -O bam -o ${out}.sort.bam ${out}.sam
rm ${out}.sam

# 4. BAM indexation
samtools index ${out}.sort.bam

# 5. Filtering: leaving only primary alignments
samtools view -F0x104 -@$cpus -hb ${out}.sort.bam > ${out}.pri.bam
samtools index ${out}.pri.bam

echo "Done! Primary alignments in ${out}.pri.bam"
