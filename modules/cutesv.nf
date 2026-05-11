process CUTESV {

    tag "cutesv_${label}"

    container params.sv_container

    publishDir "${params.outdir}/cutesv", mode: 'copy'

    input:
        val ready
        path ref
        path bam
        path bai
        val label

    output:
        path "cutesv_${label}.vcf", emit: vcf

    script:
    """
    set -euo pipefail

    case "${label}" in

        ont)
            INS_BIAS=100
            INS_RATIO=0.3
            DEL_BIAS=100
            DEL_RATIO=0.3
            ;;

        hifi)
            INS_BIAS=1000
            INS_RATIO=0.9
            DEL_BIAS=1000
            DEL_RATIO=0.5
            ;;

        clr)
            INS_BIAS=100
            INS_RATIO=0.3
            DEL_BIAS=200
            DEL_RATIO=0.5
            ;;

        *)
            echo "Unsupported label"
            exit 1
            ;;
    esac

    mkdir cutesv_tmp

    cuteSV \
        ${bam} \
        ${ref} \
        cutesv_${label}.vcf \
        cutesv_tmp \
        --max_cluster_bias_INS \${INS_BIAS} \
        --diff_ratio_merging_INS \${INS_RATIO} \
        --max_cluster_bias_DEL \${DEL_BIAS} \
        --diff_ratio_merging_DEL \${DEL_RATIO}
    """
}