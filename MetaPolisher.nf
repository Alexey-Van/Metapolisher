nextflow.enable.dsl = 2

include { CHECK        } from './modules/check_containers'
include { ALIGN_HIFI } from './modules/align'
include { ALIGN_ONT }  from './modules/align'
include { ALIGN_WGS } from './modules/align'
include { DEEPVARIANT      } from './modules/deepvariant'
include { PEPPER           } from './modules/pepper'
include { ORIGINAL_T2T     } from './modules/original_t2t'
include { AUTO_POLISH      } from './modules/automated_polishing'
include { NP2              } from './modules/nextpolish2'
include { SNIFFLES as SNIFFLES_ONT   } from './modules/sniffles'
include { SNIFFLES as SNIFFLES_HIFI  } from './modules/sniffles'
include { CUTESV as CUTESV_ONT       } from './modules/cutesv'
include { CUTESV as CUTESV_HIFI      } from './modules/cutesv'
include { FLAGGER          } from './modules/flagger'
include { MERQURY          } from './modules/merqury'
include { MEDAKA           } from './modules/medaka'
include { MERGE_READS } from './modules/merge_reads'

workflow {
    ready       = CHECK()
    draft       = Channel.fromPath(params.draft)
    draft_fai = Channel.fromPath(params.draft+".fai")

    ont_files = Channel.fromPath(
        "${params.ont}/**/*.fq.gz",
        checkIfExists: true
    )

    ont = ont_files
        .collect()
        .map { files -> tuple("ont", files) }

    ont = MERGE_READS(ont)

    
    wgs_r1_files = Channel.empty()
    wgs_r2_files = Channel.empty()

    if (params.illumina) {

        wgs_r1_files = Channel.fromPath(
            "${params.illumina}/**/*1.fq.gz",
            checkIfExists: true
        )

        wgs_r2_files = Channel.fromPath(
            "${params.illumina}/**/*2.fq.gz",
            checkIfExists: true
        )

    } else if (params.mgi) {

        wgs_r1_files = Channel.fromPath(
            "${params.mgi}/**/*1.fq.gz",
            checkIfExists: true
        )

        wgs_r2_files = Channel.fromPath(
            "${params.mgi}/**/*2.fq.gz",
            checkIfExists: true
        )

    } else {

        error "ERROR: Provide --illumina or --mgi"
    }

    wgs_r1 = wgs_r1_files
        .collect()
        .map { files -> tuple("wgs_R1", files) }

    wgs_r2 = wgs_r2_files
        .collect()
        .map { files -> tuple("wgs_R2", files) }

    wgs_r1 = MERGE_READS(wgs_r1)
    wgs_r2 = MERGE_READS(wgs_r2)

    align_illumina = ALIGN_WGS( 
        ready, 
        draft, 
        wgs_r1, 
        wgs_r2 
    )

    align_ont = ALIGN_ONT( 
        ready, 
        draft, 
        ont 
    )

    if (params.hifi) {
        
        hifi_files = Channel.fromPath(
            "${params.hifi}/**/*.fq.gz",
            checkIfExists: true
        )

        hifi = hifi_files
            .collect()
            .map { files -> tuple("hifi", files) }

        hifi = MERGE_READS(hifi)

        align_hifi = ALIGN_HIFI( 
            ready, 
            draft, 
            hifi 
        )

        deepvariant = DEEPVARIANT(
            ready,
            draft,
            draft_fai,
            align_illumina.bam,
            align_illumina.bai,
            params.hifi_bam ?: "none"
        )

        pepper = PEPPER(
            ready,
            align_ont.bam,
            align_ont.bai,
            draft,
            draft_fai
        )

        t2t = ORIGINAL_T2T(
            ready,
            deepvariant.vcf,
            pepper.vcf,
            draft,
            wgs_r1,
            wgs_r2,
            hifi
        )

        auto_polish = AUTO_POLISH(
            ready,
            draft,
            hifi,
            t2t.readmers_meryl,
            "hifi"
        )

        medaka = MEDAKA(
            ready,
            draft,
            ont
        )

        np2 = NP2(
            ready,
            draft,
            align_hifi.bam,
            wgs_r1,
            wgs_r2
        )

        sniffles_ont = SNIFFLES_ONT(
            ready,
            draft,
            align_ont.bam,
            align_ont.bai,
            "ont"
        )

        sniffles_hifi = SNIFFLES_HIFI(
            ready,
            draft,
            align_hifi.bam,
            align_hifi.bai,
            "hifi"
        )

        cutesv_ont = CUTESV_ONT(
            ready,
            draft,
            align_ont.bam,
            align_ont.bai,
            "ont"
        )

        cutesv_hifi = CUTESV_HIFI(
            ready,
            draft,
            align_hifi.bam,
            align_hifi.bai,
            "hifi"
        )

        flagger = FLAGGER(
            ready,
            draft,
            align_hifi.bam
        )

        merqury = MERQURY(
            ready,
            draft,
            wgs_r1,
            wgs_r2
        )

    } else {

        deepvariant = DEEPVARIANT(
            ready,
            draft,
            draft_fai,
            align_illumina.bam,
            align_illumina.bai,
            params.hifi_bam ?: "none"
        )

        pepper = PEPPER(
            ready,
            align_ont.bam,
            align_ont.bai,
            draft,
            draft_fai
        )

        medaka = MEDAKA(
            ready,
            draft,
            ont
        )

        np2 = NP2(
            ready,
            draft,
            align_ont.bam,
            wgs_r1,
            wgs_r2
        )

        sniffles_ont = SNIFFLES_ONT(
            ready,
            draft,
            align_ont.bam,
            align_ont.bai,
            "ont"
        )

        cutesv_ont = CUTESV_ONT(
            ready,
            draft,
            align_ont.bam,
            align_ont.bai,
            "ont"
        )

        flagger = FLAGGER(
            ready,
            draft,
            align_ont.bam
        )

        merqury = MERQURY(
            ready,
            draft,
            wgs_r1,
            wgs_r2
        )

        t2t = ORIGINAL_T2T(
            ready,
            deepvariant.vcf,
            pepper.vcf,
            draft,
            wgs_r1,
            wgs_r2,
            params.hifi ?: "none"
        )

        auto_polish = AUTO_POLISH(
            ready,
            draft,
            ont,
            t2t.readmers_meryl,
            "ont"
        )

    }

}

