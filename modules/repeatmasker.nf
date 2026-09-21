process REPEATMASKER {

    publishDir "${params.outdir}/repeatmasker", mode:'copy'

    input:
    path genome

    output:
    path "repeatmasker_output.gff3"
    path "**"

    script:
    """
    ./Repeatmasker.sh ${params.repeatmasker_path} ${genome}
    """
}