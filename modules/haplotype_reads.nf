process HAPLOTYPE_READS_HIFI {

    tag "HIFI"

    container params.mapping_container

    input:
        path hap1
        path hap2
        path classifier
        path reads, arity: '1..*'

    output:
        tuple val("hap1"), path("haplotype/hifi/hifi.hap1.fastq.gz"), emit: hap1
        tuple val("hap2"), path("haplotype/hifi/hifi.hap2.fastq.gz"), emit: hap2
        path "**"

    script:
    def read_args = reads.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    mkdir -p haplotype/hifi

    cat ${read_args} > haplotype/hifi/hifi_reads.fastq.gz

    (
        minimap2 \
            -t ${task.cpus.intdiv(2)} \
            -ax map-hifi \
            --secondary=no \
            ${hap1} \
            haplotype/hifi/hifi_reads.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/hifi/hifi_hap1.bam -
    ) &
    pid1=$!

    (
        minimap2 \
            -t ${task.cpus.intdiv(2)} \
            -ax map-hifi \
            --secondary=no \
            ${hap2} \
            haplotype/hifi/hifi_reads.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/hifi/hifi_hap2.bam -
    ) &
    pid2=$!

    status=0
    wait "$pid1" || status=$?
    wait "$pid2" || status=$?
    if [ "$status" -ne 0 ]; then
        echo "ERROR: HiFi alignment failed" >&2
        exit "$status"
    fi

    python3 ${classifier} \
        --bam1 haplotype/hifi/hifi_hap1.bam \
        --bam2 haplotype/hifi/hifi_hap2.bam \
        --reads haplotype/hifi/hifi_reads.fastq.gz \
        --out-prefix haplotype/hifi/hifi \
        --min-mapq ${params.haplotype_min_mapq} \
        --min-diff 5 \
        --gzip


    zcat haplotype/hifi/hifi.hap1.fastq.gz haplotype/hifi/hifi.shared.fastq.gz \
        | gzip -c > haplotype/hifi/hifi.hap1.fastq.gz

    zcat haplotype/hifi/hifi.hap2.fastq.gz haplotype/hifi/hifi.shared.fastq.gz \
        | gzip -c > haplotype/hifi/hifi.hap2.fastq.gz
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
        tuple val("hap1"), path("haplotype/ont/ont.hap1.fastq.gz"), emit: hap1
        tuple val("hap2"), path("haplotype/ont/ont.hap2.fastq.gz"), emit: hap2
        path "**"

    script:
    def read_args = reads.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    mkdir -p haplotype/ont

    cat ${read_args} > haplotype/ont/ont_reads.fastq.gz

    (
        minimap2 \
            -t ${task.cpus.intdiv(2)} \
            -ax map-ont \
            --secondary=no \
            ${hap1} \
            haplotype/ont/ont_reads.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/ont/ont_hap1.bam -
    ) &
    pid1=$!

    (
        minimap2 \
            -t ${task.cpus.intdiv(2)} \
            -ax map-ont \
            --secondary=no \
            ${hap2} \
            haplotype/ont/ont_reads.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/ont/ont_hap2.bam -
    ) &
    pid2=$!

    status=0
    wait "$pid1" || status=$?
    wait "$pid2" || status=$?
    if [ "$status" -ne 0 ]; then
        echo "ERROR: ONT alignment failed" >&2
        exit "$status"
    fi

    python3 ${classifier} \
        --bam1 haplotype/ont/ont_hap1.bam \
        --bam2 haplotype/ont/ont_hap2.bam \
        --reads haplotype/ont/ont_reads.fastq.gz \
        --out-prefix haplotype/ont/ont \
        --min-mapq ${params.haplotype_min_mapq} \
        --min-diff 10 \
        --gzip


    zcat haplotype/ont/ont.hap1.fastq.gz haplotype/ont/ont.shared.fastq.gz \
        | gzip -c > haplotype/ont/ont.hap1.fastq.gz

    zcat haplotype/ont/ont.hap2.fastq.gz haplotype/ont/ont.shared.fastq.gz \
        | gzip -c > haplotype/ont/ont.hap2.fastq.gz
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
        tuple val("hap1"), path("haplotype/wgs/wgs.hap1.R1.fastq.gz"), path("haplotype/wgs/wgs.hap1.R2.fastq.gz"), emit: hap1
        tuple val("hap2"), path("haplotype/wgs/wgs.hap2.R1.fastq.gz"), path("haplotype/wgs/wgs.hap2.R2.fastq.gz"), emit: hap2
        path "**"

    script:
    def r1_args = reads_r1.collect { it.toString() }.join(' ')
    def r2_args = reads_r2.collect { it.toString() }.join(' ')

    """
    set -euo pipefail

    mkdir -p haplotype/wgs

    cat ${r1_args} > haplotype/wgs/wgs_reads_R1.fastq.gz
    cat ${r2_args} > haplotype/wgs/wgs_reads_R2.fastq.gz

    bwa-mem2 index ${hap1}
    bwa-mem2 index ${hap2}

    (
        bwa-mem2 mem \
            -t ${task.cpus.intdiv(2)} \
            --secondary=no \
            ${hap1} \
            haplotype/wgs/wgs_reads_R1.fastq.gz \
            haplotype/wgs/wgs_reads_R2.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/wgs/wgs_hap1.bam -
    ) &
    pid1=$!

    (
        bwa-mem2 mem \
            -t ${task.cpus.intdiv(2)} \
            --secondary=no \
            ${hap2} \
            haplotype/wgs/wgs_reads_R1.fastq.gz \
            haplotype/wgs/wgs_reads_R2.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/wgs/wgs_hap2.bam -
    ) &
    pid2=$!

    status=0
    wait "$pid1" || status=$?
    wait "$pid2" || status=$?
    if [ "$status" -ne 0 ]; then
        echo "ERROR: WGS alignment failed" >&2
        exit "$status"
    fi

    python3 ${classifier} \
        --bam1 haplotype/wgs/wgs_hap1.bam \
        --bam2 haplotype/wgs/wgs_hap2.bam \
        --r1 haplotype/wgs/wgs_reads_R1.fastq.gz \
        --r2 haplotype/wgs/wgs_reads_R2.fastq.gz \
        --out-prefix haplotype/wgs/wgs \
        --min-mapq ${params.haplotype_min_mapq} \
        --gzip

    zcat haplotype/wgs/wgs.hap1.R1.fq.gz haplotype/wgs/wgs.shared.R1.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap1.R1.fastq.gz

    zcat haplotype/wgs/wgs.hap1.R2.fq.gz haplotype/wgs/wgs.shared.R2.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap1.R2.fastq.gz

    zcat haplotype/wgs/wgs.hap2.R1.fq.gz haplotype/wgs/wgs.shared.R1.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap2.R1.fastq.gz
    
    zcat haplotype/wgs/wgs.hap2.R2.fq.gz haplotype/wgs/wgs.shared.R2.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap2.R2.fastq.gz
    """
}