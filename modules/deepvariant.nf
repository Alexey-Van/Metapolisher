process DEEPVARIANT {

    tag "deepvariant"
    container params.deepvariant_container
    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
        val ready
        path ref
        path ref_fai
        path illumina_bam
        path illumina_bai
        val hifi_bam   // <-- теперь val, а не path

    output:
        path "deepvariant.vcf.gz", emit: vcf
        path "deepvariant.vcf.gz.tbi"
        path "**"

    script:
    """
    set -euo pipefail

    if [ "${hifi_bam}" != "none" ] && [ -f "${hifi_bam}" ]; then
        echo "[DeepVariant] Hybrid mode: Illumina + HiFi"

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

    else
        echo "[DeepVariant] Illumina-only mode"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type WGS \
            --ref ${ref} \
            --reads ${illumina_bam} \
            --output_vcf deepvariant.vcf.gz \
            --num_shards ${task.cpus}
    fi

    bcftools index deepvariant.vcf.gz
    """
}
