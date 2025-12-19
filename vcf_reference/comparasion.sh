bgzip -c data/T2T/original_fasta/truth.fixed.vcf > data/T2T/original_fasta/truth.fixed.vcf.gz
tabix -p vcf data/T2T/original_fasta/truth.fixed.vcf.gz
bcftools isec -p compare_out \
  data/T2T/original_fasta/truth.fixed.vcf.gz \
  ml_filtered_7_trheshold_5_depth_10.norm.vcf.gz
bcftools consensus -f original.fasta ml_filtered.vcf > polished.fasta