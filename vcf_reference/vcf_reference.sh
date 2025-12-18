# 1. Индексируем оригинальный FASTA
samtools faidx /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/modified_fasta/chm13_chr1_mutated.fa

# 2. Сравниваем мутированный FASTA с оригиналом
minimap2 --cs -x asm5 /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/modified_fasta/chm13_chr1_mutated.fa /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/chm13_chr1.fa  > /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/diff.paf

# 3. Вызываем варианты напрямую из PAF
k8 /home/alex/miniconda3/envs/bcf/bin/paftools.js call /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/diff.paf > /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_raw.vcf

bcftools norm \
  -f /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/chm13_chr1.fa \
  -m -any \
  /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_or.raw.vcf \
  -Oz -o /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_or.vcf.gz

bcftools index /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_or.vcf.gz

bcftools stats /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_or.vcf.gz
