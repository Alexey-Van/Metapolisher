process HAPLOTYPE_READS_HIFI {

    tag "HIFI"

    container params.mapping_container

    input:
        path hap1, stageAs: 'hap1.fasta'
        path hap2, stageAs: 'hap2.fasta'
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

    for read_file in ${read_args}; do
        case "\$read_file" in
            *.gz) gzip -cd -- "\$read_file" ;;
            *) cat -- "\$read_file" ;;
        esac
    done | gzip -c > haplotype/hifi/hifi_reads.fastq.gz

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
        | gzip -c > haplotype/hifi/hifi.hap1.merged.fastq.gz
    mv haplotype/hifi/hifi.hap1.merged.fastq.gz haplotype/hifi/hifi.hap1.fastq.gz

    zcat haplotype/hifi/hifi.hap2.fastq.gz haplotype/hifi/hifi.shared.fastq.gz \
        | gzip -c > haplotype/hifi/hifi.hap2.merged.fastq.gz
    mv haplotype/hifi/hifi.hap2.merged.fastq.gz haplotype/hifi/hifi.hap2.fastq.gz
    """
}


process HAPLOTYPE_READS_ONT {

    tag "ONT"

    container params.mapping_container

    input:
        path hap1, stageAs: 'hap1.fasta'
        path hap2, stageAs: 'hap2.fasta'
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

    for read_file in ${read_args}; do
        case "\$read_file" in
            *.gz) gzip -cd -- "\$read_file" ;;
            *) cat -- "\$read_file" ;;
        esac
    done | gzip -c > haplotype/ont/ont_reads.fastq.gz

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
        | gzip -c > haplotype/ont/ont.hap1.merged.fastq.gz
    mv haplotype/ont/ont.hap1.merged.fastq.gz haplotype/ont/ont.hap1.fastq.gz

    zcat haplotype/ont/ont.hap2.fastq.gz haplotype/ont/ont.shared.fastq.gz \
        | gzip -c > haplotype/ont/ont.hap2.merged.fastq.gz
    mv haplotype/ont/ont.hap2.merged.fastq.gz haplotype/ont/ont.hap2.fastq.gz
    """
}


process HAPLOTYPE_READS_WGS {

    tag "WGS"

    container params.mapping_container

    input:
        path hap1, stageAs: 'hap1.fasta'
        path hap2, stageAs: 'hap2.fasta'
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

    for read_file in ${r1_args}; do
        case "\$read_file" in
            *.gz) gzip -cd -- "\$read_file" ;;
            *) cat -- "\$read_file" ;;
        esac
    done | gzip -c > haplotype/wgs/wgs_reads_R1.fastq.gz
    for read_file in ${r2_args}; do
        case "\$read_file" in
            *.gz) gzip -cd -- "\$read_file" ;;
            *) cat -- "\$read_file" ;;
        esac
    done | gzip -c > haplotype/wgs/wgs_reads_R2.fastq.gz

    bwa-mem3 index ${hap1}
    bwa-mem3 index ${hap2}

    (
        bwa-mem3 mem \
            -t ${task.cpus.intdiv(2)} \
            ${hap1} \
            haplotype/wgs/wgs_reads_R1.fastq.gz \
            haplotype/wgs/wgs_reads_R2.fastq.gz \
        | samtools sort -n -@ 1 -o haplotype/wgs/wgs_hap1.bam -
    ) &
    pid1=$!

    (
        bwa-mem3 mem \
            -t ${task.cpus.intdiv(2)} \
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
        | gzip -c > haplotype/wgs/wgs.hap1.R1.merged.fastq.gz
    mv haplotype/wgs/wgs.hap1.R1.merged.fastq.gz haplotype/wgs/wgs.hap1.R1.fastq.gz

    zcat haplotype/wgs/wgs.hap1.R2.fq.gz haplotype/wgs/wgs.shared.R2.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap1.R2.merged.fastq.gz
    mv haplotype/wgs/wgs.hap1.R2.merged.fastq.gz haplotype/wgs/wgs.hap1.R2.fastq.gz

    zcat haplotype/wgs/wgs.hap2.R1.fq.gz haplotype/wgs/wgs.shared.R1.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap2.R1.merged.fastq.gz
    mv haplotype/wgs/wgs.hap2.R1.merged.fastq.gz haplotype/wgs/wgs.hap2.R1.fastq.gz
    
    zcat haplotype/wgs/wgs.hap2.R2.fq.gz haplotype/wgs/wgs.shared.R2.fq.gz \
        | gzip -c > haplotype/wgs/wgs.hap2.R2.merged.fastq.gz
    mv haplotype/wgs/wgs.hap2.R2.merged.fastq.gz haplotype/wgs/wgs.hap2.R2.fastq.gz
    """
}
