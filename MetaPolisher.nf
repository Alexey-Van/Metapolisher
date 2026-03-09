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
// BUILD ENVIRONMENT
// ----------------------

process BUILD_ENV {

    container 'docker:latest'
    publishDir "${params.outdir}/env", mode: 'copy'
    tag "docker-build"

    input:
    path dockerfile

    output:
    path "docker_build.log"

    script:
    """
    echo ">>> Building metapolisher-tools:latest..."
    docker build -t metapolisher-tools:latest -f ${dockerfile} .

    echo ">>> Pulling external images..."
    docker pull google/deepvariant:1.9.0
    docker pull kishwars/pepper_deepvariant:r0.8
    docker pull mobinasri/flagger:v1.1.0--15d859f71ec26384837dee0add731b50aac158db

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
    tuple path(reads), val(flag)

    output:
    tuple val(flag), path("*.bam")

    script:
    """
    bash /scripts/align.sh ${flag} ${assembly} ${reads}
    """
}

// ----------------------
// TOOL PROCESSES
// ----------------------

process PARLIAMENT2 {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input:
    path bam

    output:
    path "*.vcf"

    script:
    """
    bash /scripts/parliament2.sh ${bam} ${bam}.bai ${params.assembly} ${params.assembly}.fai parliament2_out
    """
}

process DEEPVARIANT {

    container 'google/deepvariant:1.9.0'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

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
        --output_vcf=deepvariant_out.vcf.gz \
        --num_shards=8
    """
}

process PEPPER {

    container 'kishwars/pepper_deepvariant:r0.8'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "${bam.simpleName}"

    input:
    path bam

    output:
    path "*.vcf"

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

    input:
    path bam

    output:
    path "*.vcf"

    script:
    """
    bash /scripts/sniffles.sh ${params.assembly} ${bam} sniffles_out.vcf 8 3
    """
}

process JASMINE {

    container 'metapolisher-tools:latest'
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    tag "jasmine"

    input:
    path vcfs

    output:
    path "*.vcf"

    script:
    """
    bash /scripts/jasmine.sh ${vcfs}
    """
}

// ----------------------
// WORKFLOW
// ----------------------

workflow {

    assembly_ch    = Channel.fromPath(params.assembly)
    truth_vcf_ch   = Channel.fromPath(params.truth_vcf)
    annotation_ch  = Channel.fromPath(params.annotation)
    dockerfile_ch  = Channel.fromPath("${baseDir}/Dockerfile")

    reads_ch = Channel
        .fromPath(params.reads_table)
        .splitCsv(header:true, sep:' ')
        .filter { it.reads_path }
        .map { row -> tuple(file(row.reads_path), row.flag) }

    BUILD_ENV(dockerfile_ch)

    align_out = ALIGN(assembly_ch, reads_ch)

    illumina_bam = align_out.filter { flag, bam -> flag == 'illumina' }
    hifi_bam     = align_out.filter { flag, bam -> flag == 'hifi' }
    ont_bam      = align_out.filter { flag, bam -> flag == 'ont' }
    clr_bam      = align_out.filter { flag, bam -> flag == 'clr' }

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

    all_sniffles_vcf = Channel
        .merge(sniffles_hifi_vcf, sniffles_ont_vcf, sniffles_clr_vcf)
        .collect()

    jasmine_vcf = JASMINE(all_sniffles_vcf)
}