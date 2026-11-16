process DEEPVARIANT {

    tag "deepvariant"
    container params.deepvariant_container
    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
        val ready
        path ref
        path ref_fai
        path illumina_bam, stageAs: 'wgs.bam'
        path illumina_bai, stageAs: 'wgs.bam.bai'
        path hifi_bam, stageAs: 'hifi.bam'
        val has_hifi

    output:
        path "deepvariant.vcf.gz", emit: vcf
        path "deepvariant.vcf.gz.tbi"
        path "**"

    script:
    """
    set -euo pipefail

    if [ "${has_hifi}" = "true" ]; then
        echo "[DeepVariant] Hybrid mode: Illumina + HiFi"

        samtools merge \
            -@ ${task.cpus} \
            merged.bam \
            hifi.bam \
            wgs.bam

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
            --reads wgs.bam \
            --output_vcf deepvariant.vcf.gz \
            --num_shards ${task.cpus}
    fi

    bcftools index deepvariant.vcf.gz
    """
}
