process SNIFFLES {

    tag "sniffles_${label}"

    container params.sv_container

    publishDir "${params.outdir}/sniffles", mode: 'copy'

    input:
        val ready
        path ref
        path bam
        path bai
        val label

    output:
        path "sniffles_${label}.vcf", emit: vcf

    script:
    """
    set -euo pipefail

    sniffles \
        --input ${bam} \
        --vcf sniffles_${label}.vcf \
        --reference ${ref} \
        --threads ${task.cpus} \
        --minsupport 5
    """
}