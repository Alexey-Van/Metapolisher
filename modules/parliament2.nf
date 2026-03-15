process PARLIAMENT2 {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode:'copy'

    input:
    path bam

    output:
    path "*.vcf"

    script:
    """
    bash /scripts/parliament2.sh ${bam} ${bam}.bai ${params.assembly} ${params.assembly}.fai parliament2_out
    """
}