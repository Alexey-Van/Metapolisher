process MERGE_READS {

    tag "${prefix}"

    container params.mapping_container

    input:
        tuple val(prefix), path(reads)

    output:
        path "${prefix}.fq.gz"

    script:
    """
    set -euo pipefail

    cat ${reads.join(' ')} > ${prefix}.fq.gz
    """
}