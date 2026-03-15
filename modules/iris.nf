process IRIS {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/iris", mode:'copy'

    input:
    path assembly
    path bam
    path vcf

    output:
    path "*"

    script:
    """
    ./iris.sh iris_input iris_output ${assembly} ${bam} ${vcf} parliament_output
    """
}