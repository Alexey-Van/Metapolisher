process CATBOOST {

    publishDir "${params.outdir}/catboost", mode:'copy'
    tag "${flag}"

    input:
    val flag
    path vcf
    path repeat_gff
    path flagger
    path merfin

    output:
    path "*.tsv"
    path "*.vcf"

    script:
    """
    python /scripts/vcf_catboost.py \
        --vcfs ${vcf} \
        --truth_vcf ${params.truth_vcf} \
        --repeat_gff ${repeat_gff} \
        --liftoff_gff ${params.liftoff_gff} \
        --low_complex ${params.low_complex} \
        --flagger ${flagger} \
        --merfin_pass_vcf ${merfin}/merfin_pass.vcf \
        --merfin_scores_bed ${merfin}/kmer_scores.bed \
        --out_vcf ${flag}.filtered.vcf \
        --threshold 0.9 \
        --out_table ${flag}_features.tsv
    """
}