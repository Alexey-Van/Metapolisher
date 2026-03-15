process FLAGGER {

    publishDir "${params.outdir}/flagger", mode:'copy'

    input:
    path genome
    path bam

    output:
    path "flagger_prediction.bed"

    script:
    """
    ./flagger.sh workdir ${genome} ${bam}
    """
}