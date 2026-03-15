process JASMINE {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode:'copy'

    input:
    path vcfs

    output:
    path "merged.vcf"

    script:
    """
    bash /scripts/jasmine.sh ${vcfs}
    """
}