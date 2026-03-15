process MERFIN {

    publishDir "${params.outdir}/merfin", mode:'copy'

    input:
    path genome
    path vcf

    output:
    path "merfin_pass.vcf"
    path "kmer_scores.bed"

    script:
    """
    ./merfin.sh ${genome} ${params.meryl_db} ${vcf}
    """
}