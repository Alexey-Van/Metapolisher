process HARD_TRAIN {

    tag "hard_train"

    container params.ml_container

    publishDir "${params.outdir}/ml", mode: 'copy'

    input:
        path dv_vcf
        path pepper_vcf
        path t2t_vcf
        path np2_vcf
        path sniffles_ont_vcf
        path sniffles_hifi_vcf
        path cutesv_ont_vcf
        path cutesv_hifi_vcf
        path truth_vcf
        path flagger_bed
        path merqury_table

    output:
        path "*"

    script:
    """
    set -euo pipefail

    python /opt/metapolisher/scripts/train_variants.py \
        --vcfs \
            ${dv_vcf} \
            ${pepper_vcf} \
            ${t2t_vcf} \
            ${np2_vcf} \
            ${sniffles_ont_vcf} \
            ${sniffles_hifi_vcf} \
            ${cutesv_ont_vcf} \
            ${cutesv_hifi_vcf} \
        --truth_vcf ${truth_vcf} \
        --flagger ${flagger_bed} \
        --merqury ${merqury_table} \
        --prefix metapolisher
    """
}