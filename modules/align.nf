process ALIGN_HIFI {

    tag "align_hifi"

    container params.mapping_container

    publishDir "${params.outdir}/align/hifi", mode: 'copy'

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

    publishDir "${params.outdir}/align/ont", mode: 'copy'

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

process ALIGN_ILLUMINA {

    tag "align_illumina"

    container params.mapping_container

    publishDir "${params.outdir}/align/illumina", mode: 'copy'

    input:
        val ready
        path draft
        path r1
        path r2

    output:
        path "illumina.sorted.bam", emit: bam
        path "illumina.sorted.bam.bai", emit: bai

    script:
    """
    set -euo pipefail

    bwa-mem2 index ${draft}

    bwa-mem2 mem \
        -t ${task.cpus} \
        ${draft} \
        ${r1} \
        ${r2} \
        > illumina.sam

    samtools sort \
        -@ ${task.cpus} \
        -o illumina.sorted.bam \
        illumina.sam

    samtools index illumina.sorted.bam
    """
}

workflow ALIGN_ALL {

    take:
        ready
        draft
        hifi
        ont
        illumina_r1
        illumina_r2

    main:

        hifi_align = ALIGN_HIFI(ready, draft, hifi)

        ont_align = ALIGN_ONT(ready, draft, ont)

        illumina_align = ALIGN_ILLUMINA(
            ready,
            draft,
            illumina_r1,
            illumina_r2
        )

    emit:
        hifi_bam     = hifi_align.bam
        hifi_bai     = hifi_align.bai

        ont_bam      = ont_align.bam
        ont_bai      = ont_align.bai

        illumina_bam = illumina_align.bam
        illumina_bai = illumina_align.bai
}