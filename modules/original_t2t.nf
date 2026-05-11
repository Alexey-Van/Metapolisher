process ORIGINAL_T2T {

    tag "original_t2t"

    container params.polishing_container

    publishDir "${params.outdir}/t2t_merfin", mode: 'copy'

    input:
        val ready
        path deepvariant_vcf
        path pepper_vcf
        path draft
        path illumina_r1
        path illumina_r2
        path hifi

    output:
        path "merfin_results*.vcf", emit: vcf
        path "*.combined_readmers.meryl", emit: readmers_meryl
        // path "merfin_results.fa", emit: fa

    script:
    """
    set -euo pipefail

    PREFIX=t2t

    THREADS=${task.cpus}
    K=21

    echo "[1] Filtering DeepVariant"

    bcftools view \
        -i 'FORMAT/GQ>=30' \
        ${deepvariant_vcf} \
        -Oz \
        -o \${PREFIX}.deepvariant.filtered.vcf.gz

    bcftools index \${PREFIX}.deepvariant.filtered.vcf.gz

    echo "[2] Filtering PEPPER"

    bcftools view \
        -i 'FORMAT/GQ>=25' \
        ${pepper_vcf} \
        -Oz \
        -o \${PREFIX}.pepper.filtered.vcf.gz

    bcftools index \${PREFIX}.pepper.filtered.vcf.gz

    echo "[3] Normalization"

    bcftools norm \
        -f ${draft} \
        -m-any \
        \${PREFIX}.deepvariant.filtered.vcf.gz \
        -Oz \
        -o \${PREFIX}.deepvariant.norm.vcf.gz

    bcftools norm \
        -f ${draft} \
        -m-any \
        \${PREFIX}.pepper.filtered.vcf.gz \
        -Oz \
        -o \${PREFIX}.pepper.norm.vcf.gz

    bcftools index \${PREFIX}.deepvariant.norm.vcf.gz
    bcftools index \${PREFIX}.pepper.norm.vcf.gz

    echo "[4] Merge"

    bcftools merge \
        \${PREFIX}.deepvariant.norm.vcf.gz \
        \${PREFIX}.pepper.norm.vcf.gz \
        -Oz \
        -o \${PREFIX}.merged.vcf.gz

    bcftools index \${PREFIX}.merged.vcf.gz

    echo "[5] Assembly kmers"

    meryl count \
        k=\${K} \
        output \${PREFIX}.seqmers.meryl \
        ${draft}

    echo "[6] Illumina kmers"

    meryl count \
        k=\${K} \
        output \${PREFIX}.illumina_readmers.meryl \
        ${illumina_r1} \
        ${illumina_r2}

    echo "[7] HiFi kmers"

    meryl count \
        k=\${K} \
        output \${PREFIX}.hifi_readmers.meryl \
        ${hifi}

    echo "[8] Combined kmers"

    meryl union-sum \
        output \${PREFIX}.combined_readmers.meryl \
        \${PREFIX}.illumina_readmers.meryl \
        \${PREFIX}.hifi_readmers.meryl

    echo "[9] Merfin"

    merfin \
        -sequence ${draft} \
        -readmers \${PREFIX}.combined_readmers.meryl \
        -seqmers \${PREFIX}.seqmers.meryl \
        -peak 106 \
        -output merfin_results \
        -threads ${task.cpus} \
        -polish \
        -vcf \${PREFIX}.merged.vcf.gz
    """
}