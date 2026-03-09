nextflow.enable.dsl = 2

// ----------------------
// PARAMETERS
// ----------------------

params.assembly      = "${baseDir}/data/assembly.fasta"
params.truth_vcf     = "${baseDir}/data/truth.vcf"
params.annotation    = "${baseDir}/data/annotation.gff"
params.reads_table   = "${baseDir}/data/reads.tsv"
params.outdir        = "${baseDir}/results"

// ----------------------
// INPUT CHANNELS
// ----------------------

Channel.fromPath(params.assembly).set { assembly_ch }
Channel.fromPath(params.truth_vcf).set { truth_vcf_ch }
Channel.fromPath(params.annotation).set { annotation_ch }
Channel.fromPath("${baseDir}/Dockerfile").set { dockerfile_ch }

Channel
    .fromPath(params.reads_table)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple( file(row.reads_path), row.flag ) }
    .set { reads_ch }

// ----------------------
// BUILD ENVIRONMENT
// ----------------------

process BUILD_ENV {
    container 'docker:latest'
    publishDir "${params.outdir}/env", mode: 'copy'
    tag "docker-build"

    input:
        path dockerfile from dockerfile_ch

    output:
        path "docker_build.log"

    script:
        """
        echo ">>> Building metapolisher-tools:latest..."
        docker build -t metapolisher-tools:latest -f ${dockerfile} .

        echo ">>> Pulling external images..."
        docker pull google/deepvariant:1.9.0
        docker pull kishwars/pepper_deepvariant:r0.8
        docker pull mobinasri/flagger:latest

        echo ">>> Done" > docker_build.log
        """
}

// ----------------------
// ALIGNMENT
// ----------------------

process ALIGN {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/alignments", mode: 'copy'
    tag "${reads.simpleName} (${flag})"

    input:
        path assembly
        tuple path(reads), val(flag) from reads_ch

    output:
        tuple val(flag), path("*.bam")

    script:
        """
        bash /scripts/align.sh ${flag} ${assembly} ${reads}
        """
}

align_out = ALIGN(assembly_ch)

// ----------------------
// SPLIT BY READ TYPE
// ----------------------

illumina_bam = align_out.filter { flag, bam -> flag == 'illumina' }
hifi_bam     = align_out.filter { flag, bam -> flag == 'hifi' }
ont_bam      = align_out.filter { flag, bam -> flag == 'ont' }
clr_bam      = align_out.filter { flag, bam -> flag == 'clr' }

// ----------------------
// TOOL PROCESSES
// ----------------------

process PARLIAMENT2 {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        bash /scripts/parliament2.sh ${bam} ${bam}.bai ${params.assembly} ${params.assembly}.fai parliament2_out
        """
}

process DEEPVARIANT {
    container 'google/deepvariant:1.9.0'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf.gz"

    script:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=${params.assembly} \
            --reads=${bam} \
            --output_vcf=deepvariant_out.vcf.gz \
            --num_shards=8
        """
}

process PEPPER {
    container 'kishwars/pepper_deepvariant:r0.8'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        run_pepper_margin_deepvariant call_variant \
            -b ${bam} \
            -f ${params.assembly} \
            -o pepper_out \
            -p pepper \
            -t 4 \
            --ont_r9_guppy5_sup
        """
}

process SNIFFLES {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        bash /scripts/sniffles.sh ${params.assembly} ${bam} sniffles_out.vcf 8 3
        """
}

process JASMINE {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "jasmine"

    input: path vcfs
    output: path "*.vcf"

    script:
        """
        bash /scripts/jasmine.sh ${vcfs}
        """
}

process FLAGGER {
    container 'mobinasri/flagger:latest'
    publishDir "${params.outdir}/flagger", mode: 'copy'
    tag "flagger"

    input: path assembly; path bam
    output: path "*"

    script:
        """
        bash /scripts/flagger.sh flagger_work ${assembly} ${bam}
        """
}

process REPEATMASKER {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/repeatmasker", mode: 'copy'
    tag "repeatmasker"

    input: path assembly
    output: path "*"

    script:
        """
        bash /scripts/Repeatmasker.sh /usr/local/bin/repeatmasker ${assembly}
        """
}

process MERFIN {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/merfin", mode: 'copy'
    tag "${vcf.simpleName}"

    input: path vcf; path truth_vcf; path assembly
    output: path "*"

    script:
        """
        bash /scripts/merfin.sh ${assembly} meryl_db ${vcf}
        """
}

process CATBOOST {
    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/catboost", mode: 'copy'
    tag "${vcf.simpleName} (${flag})"

    input: val flag; path vcf; path annotation
    output: path "*.tsv"

    script:
        """
        python /scripts/vcf_catboost.py \
            --vcfs ${vcf} \
            --truth_vcf ${params.truth_vcf} \
            --repeat_gff /results/repeatmasker/repeatmasker_output.gff3 \
            --merfin_pass_vcf /results/merfin/merfin_pass.vcf \
            --liftoff_gff liftoff.gff3 \
            --low_complex low_complexity.bed \
            --flagger /results/flagger/flagger_prediction.bed \
            --out_vcf filtered_output.vcf \
            --threshold 0.9 \
            --out_table variant_features.tsv
        """
}

// ----------------------
// WORKFLOW LOGIC
// ----------------------

workflow {
    BUILD_ENV(dockerfile_ch)

    illumina_bam_only = illumina_bam.map{ flag, bam -> bam }.into { illum1; illum2 }
    parliament_vcf = PARLIAMENT2(illum1)
    deepvariant_illumina_vcf = DEEPVARIANT(illum2)

    hifi_bam_only = hifi_bam.map{ flag, bam -> bam }.into { h1; h2 }
    sniffles_hifi_vcf = SNIFFLES(h1)
    deepvariant_hifi_vcf = DEEPVARIANT(h2)

    ont_bam_only = ont_bam.map{ flag, bam -> bam }.into { o1; o2 }
    sniffles_ont_vcf = SNIFFLES(o1)
    pepper_vcf = PEPPER(o2)

    clr_bam_only = clr_bam.map{ flag, bam -> bam }
    sniffles_clr_vcf = SNIFFLES(clr_bam_only)

    all_sniffles_vcf = Channel.merge(sniffles_hifi_vcf, sniffles_ont_vcf, sniffles_clr_vcf).collect().set { sniffles_list }
    jasmine_vcf = JASMINE(sniffles_list)

    hifi_bam_for_flagger = hifi_bam.map{ flag, bam -> bam }.first()
    flagger_out = FLAGGER(assembly_ch, hifi_bam_for_flagger)

    repeatmasker_out = REPEATMASKER(assembly_ch)

    all_vcf = Channel.merge(parliament_vcf, deepvariant_illumina_vcf, deepvariant_hifi_vcf, pepper_vcf, jasmine_vcf)
    vcf_for_merfin = all_vcf.filter { vcf -> vcf.name != file(params.truth_vcf).name }
}