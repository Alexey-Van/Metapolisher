process DEEPVARIANT {

    tag "deepvariant"

    container params.deepvariant_container

    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
        val ready
        path ref
        path ref_fai
        path hifi_bam
        path illumina_bam

    output:
        path "deepvariant.vcf.gz", emit: vcf
        path "deepvariant.vcf.gz.tbi"

    script:
    """
    set -euo pipefail

    samtools merge \
        -@ ${task.cpus} \
        merged.bam \
        ${hifi_bam} \
        ${illumina_bam}

    samtools index merged.bam

    /opt/deepvariant/bin/run_deepvariant \
        --model_type HYBRID_PACBIO_ILLUMINA \
        --ref ${ref} \
        --reads merged.bam \
        --output_vcf deepvariant.vcf.gz \
        --num_shards ${task.cpus}

    bcftools index deepvariant.vcf.gz
    """
}