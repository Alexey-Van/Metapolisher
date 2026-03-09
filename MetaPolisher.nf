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
// GLOBAL CONTAINER
// ----------------------

process.container = 'metapolisher-tools:latest'

// ----------------------
// INPUT CHANNELS
// ----------------------

Channel.fromPath(params.assembly).set { assembly_ch }
Channel.fromPath(params.truth_vcf).set { truth_vcf_ch }
Channel.fromPath(params.annotation).set { annotation_ch }

Channel
    .fromPath(params.reads_table)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple( file(row.reads_path), row.flag ) }
    .set { reads_ch }

// ----------------------
// ALIGNMENT
// ----------------------

process ALIGN {
    publishDir "${params.outdir}/alignments", mode: 'copy'
    tag "${reads.simpleName} (${flag})"

    input:
        path assembly
        tuple path(reads), val(flag) from reads_ch

    output:
        tuple val(flag), path("*.bam")

    script:
        """
        align.sh ${flag} ${assembly} ${reads}
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
process BUILD_ENV {
    publishDir "${params.outdir}/env", mode: 'copy'
    tag "docker-build"

    input:
        path "Dockerfile"

    output:
        path "docker_build.log"

    script:
        """
        echo ">>> Building metapolisher-tools:latest..."
        docker build -t metapolisher-tools:latest -f Dockerfile .

        echo ">>> Pulling external images..."
        docker pull google/deepvariant:1.9.0
        docker pull kishwars/pepper_deepvariant:r0.8
        docker pull mobinasri/flagger:v1.1.0--15d859f71ec26384837dee0add731b50aac158db

        echo ">>> Done" > docker_build.log
        """
}

process PARLIAMENT2 {
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        parliament2.sh ${bam} ${bam}.bai ${params.assembly} ${params.assembly}.fai parliament2_out
        """
}

process DEEPVARIANT {
    container 'google/deepvariant:1.9.0'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=${params.assembly} \
            --reads=${bam} \
            --output_vcf=deepvariant_out.vcf \
            --output_gvcf=deepvariant_out.g.vcf.gz \
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
        run_pepper_deepvariant.sh \
            --bam ${bam} \
            --fasta ${params.assembly} \
            --output_dir pepper_out \
            --output_prefix pepper
        """
}

process SNIFFLES {
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input: path bam
    output: path "*.vcf"

    script:
        """
        sniffles.sh ${params.assembly} ${bam} sniffles_out.vcf 8 3
        """
}

process JASMINE {
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "jasmine"

    input: path vcfs
    output: path "*.vcf"

    script:
        """
        jasmine.sh ${vcfs}
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
        flagger.sh flagger_work ${assembly} ${bam}
        """
}

process REPEATMASKER {
    publishDir "${params.outdir}/repeatmasker", mode: 'copy'
    tag "repeatmasker"

    input: path assembly
    output: path "*"

    script:
        """
        Repeatmasker.sh /usr/local/bin/repeatmasker ${assembly}
        """
}

process MERFIN {
    publishDir "${params.outdir}/merfin", mode: 'copy'
    tag "${vcf.simpleName}"

    input: path vcf; path truth_vcf; path assembly
    output: path "*"

    script:
        """
        merfin.sh ${assembly} meryl_db ${vcf}
        """
}

process CATBOOST {
    publishDir "${params.outdir}/catboost", mode: 'copy'
    tag "${vcf.simpleName} (${flag})"

    input: val flag; path vcf; path annotation
    output: path "*.tsv"

    script:
        """
        python vcf_catboost.py \
            --vcfs ${vcf} \
            --truth_vcf ${params.truth_vcf} \
            --repeat_gff repeatmasker_output.gff3 \
            --merfin_pass_vcf merfin_pass.vcf \
            --liftoff_gff liftoff.gff3 \
            --low_complex low_complexity.bed \
            --flagger flagger_prediction.bed \
            --out_vcf filtered_output.vcf \
            --threshold 0.9 \
            --out_table variant_features.tsv
        """
}

// ----------------------
// WORKFLOW LOGIC
// ----------------------
BUILD_ENV(Dockerfile)

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
merfin_out = MERFIN(vcf_for_merfin, truth_vcf_ch, assembly_ch)

snv_vcf = Channel.merge(deepvariant_illumina_vcf, deepvariant_hifi_vcf, pepper_vcf).map { vcf -> tuple('snv', vcf) }
sv_vcf  = Channel.merge(parliament_vcf, jasmine_vcf).map { vcf -> tuple('sv', vcf) }

all_for_catboost = Channel.merge(snv_vcf, sv_vcf)
CATBOOST(all_for_catboost.map{ flag, vcf -> tuple(flag, vcf) }, annotation_ch)

workflow.onComplete {
    println "Pipeline finished. Output directory: ${params.outdir}"
}
