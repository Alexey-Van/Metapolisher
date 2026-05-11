# 1. Вырезаем 1 Мб из BAM
samtools view -b aln.bam chr22:1-1000000 > aln.test.bam
samtools index aln.test.bam

# 2. Делаем FASTQ
samtools fastq -1 R1.test.fastq -2 R2.test.fastq aln.test.bam

# 3. Вырезаем тот же регион из VCF
bcftools view -r chr22:1-1000000 variants.vcf.gz -Oz -o variants.test.vcf.gz

# 4. Вырезаем регион из сборки
samtools faidx draft.fasta chr22:1-1000000 > draft.test.fasta
