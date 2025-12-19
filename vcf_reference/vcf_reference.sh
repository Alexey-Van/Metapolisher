samtools faidx /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/modified_fasta/chm13_chr1_mutated.fa

minimap2 --cs -x asm5 /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/modified_fasta/chm13_chr1_mutated.fa /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/chm13_chr1.fa  > /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/diff.paf

k8 /home/alex/miniconda3/envs/bcf/bin/paftools.js call /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/diff.paf > /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/true_raw.vcf

awk 'BEGIN{
  OFS="\t";
  print "##fileformat=VCFv4.2";
  print "##contig=<ID=chr1,length=248255384>";
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER";
}
$1!="R" {
  chrom=$2;
  pos=$3+1;
  ref=toupper($7);
  alt=toupper($8);
  qual=$6;              
  filter="PASS";
  print chrom, pos, ".", ref, alt, qual, filter, info, format, sample;
}' data/T2T/original_fasta/true_raw.vcf > data/T2T/original_fasta/truth.vcf

awk -v REFFA="data/T2T/modified_fasta/chm13_chr1_mutated.fa" '
BEGIN{
  OFS="\t";
  print "##fileformat=VCFv4.2";
  print "##contig=<ID=chr1,length=248255384>";
  print "##INFO=<ID=ML_PROB,Number=1,Type=Float,Description=\"CatBoost TRUE_VARIANT probability\">";
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}
function anchor_base(chrom, pos,   cmd, seq) {
  cmd = "samtools faidx " REFFA " " chrom ":" pos "-" pos " | grep -v \"^>\" | tr -d \"\\n\"";
  seq=""; cmd | getline seq; close(cmd);
  seq=toupper(seq);
  if (seq=="" || seq ~ /[^ACGTN]/) {
    printf("WARN: empty/non-ACGTN anchor at %s:%d, using N\n", chrom, pos) > "/dev/stderr";
    seq="N";
  }
  return substr(seq,1,1);
}
/^#/ { next }

{
  chrom=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5); qual=$6; filt=$7;

  if (ref=="-" && alt!="-") {
    newpos = pos - 1;
    anc = anchor_base(chrom, newpos);
    ref = anc;
    alt = anc alt;
    pos = newpos;
  } else if (alt=="-" && ref!="-") {
    newpos = pos - 1;
    anc = anchor_base(chrom, newpos);
    ref = anc ref;
    alt = anc;
    pos = newpos;
  } else if (ref=="-" && alt=="-") {
    printf("WARN: both REF and ALT are -, skipping %s:%d\n", chrom, pos) > "/dev/stderr";
    next;
  }

  info="from_paf";
  print chrom, pos, id, ref, alt, qual, filt, info;
}' data/T2T/original_fasta/truth.vcf > data/T2T/original_fasta/truth.fixed.vcf



bcftools norm -f /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/modified_fasta/chm13_chr1_mutated.fa --check-ref s \
    -Oz -o data/T2T/original_fasta/truth.norm.vcf.gz \
    data/T2T/original_fasta/truth.fixed.vcf

bcftools stats /mnt/c/Users/1/Desktop/Project_MetaPolisher/data/T2T/original_fasta/truth.norm.vcf.gz
