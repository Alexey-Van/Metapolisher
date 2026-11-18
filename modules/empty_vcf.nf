process EMPTY_VCF {

    tag "empty_vcf"

    container params.polishing_container

    input:
        path source_vcf

    output:
        path "empty.vcf", emit: vcf

    script:
    """
    set -euo pipefail

    bcftools view -h ${source_vcf} \
        | sed '/^#CHROM/ s/\tFORMAT.*\$//' \
        > empty.vcf
    """
}
