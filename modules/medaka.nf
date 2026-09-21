process MEDAKA {

    tag "medaka"

    container params.polishing_container

    publishDir "${params.outdir}/medaka", mode: 'copy'

    input:
        val ready
        path draft
        path reads

    output:
        path "medaka_output/**", emit: medaka_results

    script:
    """
    set -euo pipefail

    THREADS=${task.cpus}

    echo "[MEDAKA] Assembly: ${draft}"
    echo "[MEDAKA] ONT reads: ${reads}"
    echo "[MEDAKA] Threads: \${THREADS}"

    mkdir -p medaka_output

    medaka_consensus \
        -i ${reads} \
        -d ${draft} \
        -o medaka_output \
        -t \${THREADS}

    echo "[MEDAKA] Finished"
    """
}