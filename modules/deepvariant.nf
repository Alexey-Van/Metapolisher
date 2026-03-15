process DEEPVARIANT {

    container 'google/deepvariant:1.9.0'
    publishDir "${params.outdir}/vcfs", mode:'copy'

    input:
    path bam

    output:
    path "*.vcf.gz"

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${params.assembly} \
        --reads=${bam} \
        --output_vcf=deepvariant.vcf.gz \
        --num_shards=8
    """
}