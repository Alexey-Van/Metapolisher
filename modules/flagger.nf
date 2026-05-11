process FLAGGER {

    tag "flagger"

    container params.flagger_container

    publishDir "${params.outdir}/flagger", mode: 'copy'

    input:
        val ready
        path fasta
        path bam

    output:
        path "hmm_flagger_outputs", emit: bed

    script:
    """
    set -euo pipefail

    # 1. Индексы
    samtools faidx "${fasta}"
    samtools index "${bam}"

    # 2. Whole genome BED
    awk '{print \$1"\\t0\\t"\$2}' "${fasta}.fai" > whole_genome.bed

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
    """
}
