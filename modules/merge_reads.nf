process MERGE_READS {

    tag "${prefix}"

    container params.mapping_container

    input:
        tuple val(prefix), path(reads, arity: '1..*')

    output:
        path "${prefix}.fq.gz"

    script:
    """
    set -euo pipefail

    for read_file in ${reads.collect { it.toString() }.join(' ')}; do
        case "\$read_file" in
            *.gz) gzip -cd -- "\$read_file" ;;
            *) cat -- "\$read_file" ;;
        esac
    done | gzip -c > ${prefix}.fq.gz
    """
}
