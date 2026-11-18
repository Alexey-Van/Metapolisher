process FLAGGER {

    tag "flagger"

    container params.flagger_container

    publishDir "${params.outdir}/flagger", mode: 'copy'

    input:
        val ready
        path fasta
        path reference_indexes, arity: '1..2'
        path bam

    output:
        path "flagger.bed", emit: bed
        path "hmm_flagger_outputs", emit: results

    script:
    """
    set -euo pipefail

    # 1. Индексы
    samtools faidx "${fasta}"
    samtools index "${bam}"

    # 2. Whole genome BED
    awk '{print \$1"\t0\t"\$2}' "${fasta}.fai" > whole_genome.bed

    # 3. JSON (важно: экранируем \$)
    cat > annotations_path.json <<EOF
{
  "whole_genome": "\$(pwd)/whole_genome.bed"
}
EOF

    # 4. bam2cov
    bam2cov \
        --bam "${bam}" \
        --output coverage_file.cov.gz \
        --annotationJson annotations_path.json \
        --threads ${task.cpus} \
        --baselineAnnotation whole_genome

    # 5. hmm_flagger
    mkdir -p hmm_flagger_outputs

    hmm_flagger \
        --input coverage_file.cov.gz \
        --outputDir hmm_flagger_outputs \
        --iterations 50 \
        --threads ${task.cpus}

    # Use the final prediction track as the ML annotation. Older Flagger
    # builds use varying basenames, so prefer a final BED and fall back to
    # the largest BED in the output directory.
    candidate=\$(find hmm_flagger_outputs -maxdepth 1 -type f -name '*.bed' \
        | grep -i 'final' | sort | head -n 1 || true)
    if [ -z "\${candidate}" ]; then
        candidate=\$(find hmm_flagger_outputs -maxdepth 1 -type f -name '*.bed' \
            -printf '%s\t%p\n' | sort -nr | head -n 1 | cut -f2-)
    fi
    if [ -z "\${candidate}" ]; then
        echo "ERROR: Flagger did not produce a BED track" >&2
        exit 1
    fi
    cp "\${candidate}" flagger.bed
    """
}
