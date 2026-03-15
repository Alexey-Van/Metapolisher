process ALIGN {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/alignments", mode:'copy'

    input:
    path assembly
    tuple path(reads), val(flag)

    output:
    tuple val(flag), path("*.bam")

    script:
    """
    micromamba activate mapping 
    bash /scripts/align.sh ${flag} ${assembly} ${reads}
    """
}