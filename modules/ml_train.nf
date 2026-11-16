process ML_TRAIN {

    tag "ml_train_${prefix}"

    container params.ml_container

    publishDir "${params.outdir}/ml_train/${prefix}", mode: 'copy'

    input:
        val ready
        path script
        path vcfs, stageAs: 'vcfs/*', arity: '1..*'
        path truth_vcf
        path flagger_bed
        path merqury_bed
        val prefix

    output:
        path "model_${prefix}.cbm", emit: model
        path "${prefix}.variants.tsv", emit: table
        path "${prefix}_hyperparams.tsv", emit: hyperparameters

    script:
    def vcf_args = vcfs.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    python3 ${script} \
        --mode train \
        --vcfs ${vcf_args} \
        --truth_vcf ${truth_vcf} \
        --flagger ${flagger_bed} \
        --merqury ${merqury_bed} \
        --prefix ${prefix}
    """
}
