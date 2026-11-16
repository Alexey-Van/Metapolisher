process MEDAKA {

    tag "medaka"

    container params.polishing_container

    publishDir "${params.outdir}/medaka", mode: 'copy'

    input:
        val ready
        path draft
        path reads, arity: '1..*'

    output:
        path "medaka_output/**", emit: medaka_results

    script:
    def read_args = reads.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    THREADS=${task.cpus}

    echo "[MEDAKA] Assembly: ${draft}"
    echo "[MEDAKA] ONT reads: ${reads}"
    echo "[MEDAKA] Threads: \${THREADS}"

    mkdir -p medaka_output

    for reads_file in ${read_args}; do
        case "\$reads_file" in
            *.gz) gzip -cd -- "\$reads_file" ;;
            *) cat -- "\$reads_file" ;;
        esac
    done > ont_reads.fastq

    medaka_consensus \
        -i ont_reads.fastq \
        -d ${draft} \
        -o medaka_output \
        -t \${THREADS}

    echo "[MEDAKA] Finished"
    """
}
