process ALIGN_HIFI {

    tag "align_hifi"

    container params.mapping_container

    publishDir "${params.outdir}/align/hifi/${draft.baseName}", mode: 'copy'

    input:
        val ready
        path draft
        path hifi

    output:
        path "hifi.sorted.bam", emit: bam
        path "hifi.sorted.bam.bai", emit: bai

    script:
    """
    set -euo pipefail

    meryl count k=15 output meryl_db ${draft}

    meryl print greater-than 100 meryl_db > repetitive_k15.txt

    winnowmap \
        --MD \
        -W repetitive_k15.txt \
        -t ${task.cpus} \
        -ax map-pb \
        ${draft} \
        ${hifi} \
        > hifi.sam

    samtools sort \
        -@ ${task.cpus} \
        -o hifi.sorted.bam \
        hifi.sam

    samtools index hifi.sorted.bam
    """
}


process ALIGN_ONT {

    tag "align_ont"

    container params.mapping_container

    publishDir "${params.outdir}/align/ont/${draft.baseName}", mode: 'copy'

    input:
        val ready
        path draft
        path ont

    output:
        path "ont.sorted.bam", emit: bam
        path "ont.sorted.bam.bai", emit: bai

    script:
    """
    set -euo pipefail

    meryl count k=15 ${draft} output meryl_db

    meryl print greater-than 100 meryl_db > repetitive_k15.txt

    winnowmap \
        --MD \
        -W repetitive_k15.txt \
        -t ${task.cpus} \
        -ax map-ont \
        ${draft} \
        ${ont} \
        > ont.sam

    samtools view ont.sam -H > ont.filtered.sam

    samtools view ont.sam | awk 'length(\$6) < 64000' >> ont.filtered.sam

    samtools sort \
        -@ ${task.cpus} \
        -o ont.sorted.bam \
        ont.filtered.sam

    samtools index ont.sorted.bam
    """
}


process ALIGN_WGS {

    tag "align_wgs"

    container params.mapping_container

    publishDir "${params.outdir}/align/wgs/${draft.baseName}", mode: 'copy'

    input:
        val ready
        path draft
        path r1
        path r2

    output:
        path "wgs.sorted.bam", emit: bam
        path "wgs.sorted.bam.bai", emit: bai

    script:
    """
    set -euo pipefail

    bwa-mem2 index ${draft}

    bwa-mem2 mem \
        -t ${task.cpus} \
        ${draft} \
        ${r1} \
        ${r2} \
        > wgs.sam

    samtools sort \
        -@ ${task.cpus} \
        -o wgs.sorted.bam \
        wgs.sam

    samtools index wgs.sorted.bam
    """
}