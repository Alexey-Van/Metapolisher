#!/bin/bash

SCAFF="$1"
READS="$2"
THREADS=8
PARTS=5

mkdir -p split_reads paf_parts

# 1. Split reads
seqkit split -p $PARTS -O split_reads $READS

# 2. Run minimap2 on each chunk
for f in split_reads/*.fasta; do
    base=$(basename $f .fasta)
    minimap2 -x ava-ont -t $THREADS $SCAFF $f > paf_parts/$base.paf
done

# 3. Merge PAF
cat paf_parts/*.paf > merged.paf

echo "PAF готов: merged.paf"
