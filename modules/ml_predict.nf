process ML_PREDICT {

    tag "ml_predict_${prefix}"

    container params.ml_container

    publishDir "${params.outdir}/ml/${prefix}", mode: 'copy'

    input:
        val ready
        path script
        path model
        path vcfs, stageAs: 'vcfs/*', arity: '1..*'
        path flagger_bed
        path merqury_bed
        val prefix

    output:
        path "${prefix}.ml.vcf", emit: vcf
        path "${prefix}.ml.variants.tsv", emit: table

    script:
    def vcf_args = vcfs.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    python3 ${script} \
        --mode predict \
        --vcfs ${vcf_args} \
        --flagger ${flagger_bed} \
        --merqury ${merqury_bed} \
        --ml_model ${model} \
        --threshold ${params.ml_threshold} \
        --prefix ${prefix} \
        --out_vcf ${prefix}.ml.vcf \
        --out_table ${prefix}.ml.variants.tsv
    """
}
