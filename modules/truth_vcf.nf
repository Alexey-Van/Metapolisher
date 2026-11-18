process TRUTH_VCF {

    tag "truth_vcf_${prefix}"

    container params.truth_container

    publishDir "${params.outdir}/truth_vcf/${prefix}", mode: 'copy'

    input:
        val ready
        path truth_script
        path ref
        path ref_indexes, arity: '1..2'
        path query
        val prefix

    output:
        path "truth.vcf.gz", emit: vcf
        path "truth.vcf.gz.csi"

    script:
    """
    set -euo pipefail

    bash "${truth_script}" "${prefix}" "${ref}" "${query}" "${task.cpus}"
    """
}
