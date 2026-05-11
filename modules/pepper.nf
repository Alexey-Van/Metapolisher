process PEPPER {

    tag "pepper"

    container params.pepper_container

    publishDir "${params.outdir}/pepper", mode: 'copy'

    input:
        val ready
        path bam
        path bai
        path ref
        path fai

    output:
        path "pepper_out/*.vcf.gz", emit: vcf

    script:
    """
    set -euo pipefail

    run_pepper_margin_deepvariant call_variant \
        -b ${bam} \
        -f ${ref} \
        -o pepper_out \
        -p pepper \
        -t ${task.cpus} \
        --ont_r9_guppy5_sup
    """
}