process PEPPER {

    container 'kishwars/pepper_deepvariant:r0.8'
    publishDir "${params.outdir}/vcfs", mode:'copy'

    input:
    path bam

    output:
    path "*.vcf"

    script:
    """
    run_pepper_margin_deepvariant call_variant \
        -b ${bam} \
        -f ${params.assembly} \
        -o pepper_out \
        -p pepper \
        -t 4 \
        --ont_r9_guppy5_sup
    """
}