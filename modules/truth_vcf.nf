process TRUTH_VCF {

    tag "truth_vcf_${prefix}"

    container params.truth_container

    publishDir "${params.outdir}/truth_vcf/${prefix}", mode: 'copy'

    input:
        val ready
        path truth_script
        path ref, stageAs: 'ref.fa'
        path query, stageAs: 'query.fa'
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
