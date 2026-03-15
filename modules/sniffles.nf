process SNIFFLES {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode:'copy'

    input:
    path bam

    output:
    path "*.vcf"

    script:
    """
    bash /scripts/sniffles.sh ${params.assembly} ${bam} sniffles.vcf 8 3
    """
}