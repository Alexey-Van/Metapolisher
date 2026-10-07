process HAPLOTYPE_READS_HIFI {

    tag "HIFI"

    container params.mapping_container

    input:
        path hap1
        path hap2
        path classifier
        path reads, arity: '1..*'

    output:
        tuple val("hap1"), path("hap1.fastq.gz"), emit: hap1
        tuple val("hap2"), path("hap2.fastq.gz"), emit: hap2
        path "unassigned.fastq.gz", emit: unassigned
        path "summary.tsv", emit: summary

    script:
    def read_args = reads.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    cat ${read_args} > reads.fastq.gz

    minimap2 \
        -t ${task.cpus} \
        -ax map-hifi \
        --secondary=no \
        ${hap1} \
        reads.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap1.bam -

    minimap2 \
        -t ${task.cpus} \
        -ax map-hifi \
        --secondary=no \
        ${hap2} \
        reads.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap2.bam -

    python3 ${classifier} \
        --hap1-bam hap1.bam \
        --hap2-bam hap2.bam \
        --reads reads.fastq.gz \
        --outdir assignment \
        --delta ${params.haplotype_delta} \
        --min-mapq ${params.haplotype_min_mapq}

    zcat assignment/hap1.fastq.gz assignment/shared.fastq.gz \
        | gzip -c > hap1.fastq.gz

    zcat assignment/hap2.fastq.gz assignment/shared.fastq.gz \
        | gzip -c > hap2.fastq.gz

    cp assignment/unassigned.fastq.gz unassigned.fastq.gz
    cp assignment/summary.tsv summary.tsv
    """
}


process HAPLOTYPE_READS_ONT {

    tag "ONT"

    container params.mapping_container

    input:
        path hap1
        path hap2
        path classifier
        path reads, arity: '1..*'

    output:
        tuple val("hap1"), path("hap1.fastq.gz"), emit: hap1
        tuple val("hap2"), path("hap2.fastq.gz"), emit: hap2
        path "unassigned.fastq.gz", emit: unassigned
        path "summary.tsv", emit: summary

    script:
    def read_args = reads.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    cat ${read_args} > reads.fastq.gz

    minimap2 \
        -t ${task.cpus} \
        -ax map-ont \
        --secondary=no \
        ${hap1} \
        reads.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap1.bam -

    minimap2 \
        -t ${task.cpus} \
        -ax map-ont \
        --secondary=no \
        ${hap2} \
        reads.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap2.bam -

    python3 ${classifier} \
        --hap1-bam hap1.bam \
        --hap2-bam hap2.bam \
        --reads reads.fastq.gz \
        --outdir assignment \
        --delta ${params.haplotype_delta} \
        --min-mapq ${params.haplotype_min_mapq}

    zcat assignment/hap1.fastq.gz assignment/shared.fastq.gz \
        | gzip -c > hap1.fastq.gz

    zcat assignment/hap2.fastq.gz assignment/shared.fastq.gz \
        | gzip -c > hap2.fastq.gz

    cp assignment/unassigned.fastq.gz unassigned.fastq.gz
    cp assignment/summary.tsv summary.tsv
    """
}


process HAPLOTYPE_READS_WGS {

    tag "WGS"

    container params.mapping_container

    input:
        path hap1
        path hap2
        path classifier
        path reads_r1, arity: '1..*'
        path reads_r2, arity: '1..*'

    output:
        tuple val("hap1"), path("hap1_R1.fastq.gz"), path("hap1_R2.fastq.gz"), emit: hap1
        tuple val("hap2"), path("hap2_R1.fastq.gz"), path("hap2_R2.fastq.gz"), emit: hap2
        path "unassigned_R1.fastq.gz", emit: unassigned_r1
        path "unassigned_R2.fastq.gz", emit: unassigned_r2
        path "summary.tsv", emit: summary

    script:
    def r1_args = reads_r1.collect { it.toString() }.join(' ')
    def r2_args = reads_r2.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    cat ${r1_args} > reads_R1.fastq.gz
    cat ${r2_args} > reads_R2.fastq.gz

    minimap2 \
        -t ${task.cpus} \
        -ax sr \
        --secondary=no \
        ${hap1} \
        reads_R1.fastq.gz \
        reads_R2.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap1.bam -

    minimap2 \
        -t ${task.cpus} \
        -ax sr \
        --secondary=no \
        ${hap2} \
        reads_R1.fastq.gz \
        reads_R2.fastq.gz \
        | samtools sort -n -@ ${task.cpus} -o hap2.bam -

    python3 ${classifier} \
        --hap1-bam hap1.bam \
        --hap2-bam hap2.bam \
        --reads-r1 reads_R1.fastq.gz \
        --reads-r2 reads_R2.fastq.gz \
        --outdir assignment \
        --paired \
        --delta ${params.haplotype_delta} \
        --min-mapq ${params.haplotype_min_mapq}

    zcat assignment/hap1_R1.fastq.gz assignment/shared_R1.fastq.gz \
        | gzip -c > hap1_R1.fastq.gz

    zcat assignment/hap1_R2.fastq.gz assignment/shared_R2.fastq.gz \
        | gzip -c > hap1_R2.fastq.gz

    zcat assignment/hap2_R1.fastq.gz assignment/shared_R1.fastq.gz \
        | gzip -c > hap2_R1.fastq.gz

    zcat assignment/hap2_R2.fastq.gz assignment/shared_R2.fastq.gz \
        | gzip -c > hap2_R2.fastq.gz

    cp assignment/unassigned_R1.fastq.gz unassigned_R1.fastq.gz
    cp assignment/unassigned_R2.fastq.gz unassigned_R2.fastq.gz
    cp assignment/summary.tsv summary.tsv
    """
}