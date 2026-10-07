nextflow.enable.dsl = 2

include { CHECK        } from './modules/check_containers'
include { ALIGN_HIFI  } from './modules/align'
include { ALIGN_ONT   } from './modules/align'
include { ALIGN_WGS   } from './modules/align'
include { DEEPVARIANT } from './modules/deepvariant'
include { PEPPER      } from './modules/pepper'
include { ORIGINAL_T2T } from './modules/original_t2t'
include { AUTO_POLISH } from './modules/automated_polishing'
include { NP2         } from './modules/nextpolish2'
include { SNIFFLES as SNIFFLES_ONT  } from './modules/sniffles'
include { SNIFFLES as SNIFFLES_HIFI } from './modules/sniffles'
include { CUTESV as CUTESV_ONT  } from './modules/cutesv'
include { CUTESV as CUTESV_HIFI } from './modules/cutesv'
include { FLAGGER      } from './modules/flagger'
include { MERQURY     } from './modules/merqury'
include { MEDAKA      } from './modules/medaka'
include { MERGE_READS } from './modules/merge_reads'
include { HAPLOTYPE_READS_HIFI } from './modules/haplotype_reads'
include { HAPLOTYPE_READS_ONT  } from './modules/haplotype_reads'
include { HAPLOTYPE_READS_WGS  } from './modules/haplotype_reads'


workflow {

    ready = CHECK()


    // -----------------------------------------------------------------------
    // WGS reads
    // -----------------------------------------------------------------------

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


    // Сохраняем исходные списки для haplotype classifier.
    // MERGE_READS продолжает использоваться для обычной ветки.

    wgs_r1_raw = wgs_r1_files.collect()
    wgs_r2_raw = wgs_r2_files.collect()

    wgs_r1 = wgs_r1_files.collect().map { files ->
        tuple("wgs_R1", files)
    }

    wgs_r2 = wgs_r2_files.collect().map { files ->
        tuple("wgs_R2", files)
    }

    wgs_r1 = MERGE_READS(wgs_r1)
    wgs_r2 = MERGE_READS(wgs_r2)


    // -----------------------------------------------------------------------
    // DIPLOID
    // -----------------------------------------------------------------------

    if (params.ploidy == "diploid") {

        if (!params.hap1 || !params.hap2) {
            error "ERROR: ploidy=diploid requires --hap1 and --hap2"
        }


        // -------------------------------------------------------------------
        // Haplotype assemblies
        // -------------------------------------------------------------------

        hap1 = Channel.value(
            file(params.hap1, checkIfExists: true)
        )

        hap2 = Channel.value(
            file(params.hap2, checkIfExists: true)
        )

        hap1_fai = Channel.value(
            file("${params.hap1}.fai", checkIfExists: true)
        )

        hap2_fai = Channel.value(
            file("${params.hap2}.fai", checkIfExists: true)
        )

        classifier = Channel.value(
            file(
                "${baseDir}/scripts/assign_haplotypes.py",
                checkIfExists: true
            )
        )


        // -------------------------------------------------------------------
        // Haplotype read assignment
        // -------------------------------------------------------------------

        if (params.hifi) {

            hifi_reads = Channel
                .fromPath(params.hifi, checkIfExists: true)
                .collect()

            HAPLOTYPE_READS_HIFI(
                hap1,
                hap2,
                classifier,
                hifi_reads
            )

            hifi_hap1 = HAPLOTYPE_READS_HIFI.out.hap1.map { it[1] }
            hifi_hap2 = HAPLOTYPE_READS_HIFI.out.hap2.map { it[1] }

        }


        if (params.ont) {

            ont_reads = Channel
                .fromPath(params.ont, checkIfExists: true)
                .collect()

            HAPLOTYPE_READS_ONT(
                hap1,
                hap2,
                classifier,
                ont_reads
            )

            ont_hap1 = HAPLOTYPE_READS_ONT.out.hap1.map { it[1] }
            ont_hap2 = HAPLOTYPE_READS_ONT.out.hap2.map { it[1] }
        }


        HAPLOTYPE_READS_WGS(
            hap1,
            hap2,
            classifier,
            wgs_r1_raw,
            wgs_r2_raw
        )

        wgs_hap1 = HAPLOTYPE_READS_WGS.out.hap1.map { item ->
            [
                name: item[0],
                r1: item[1],
                r2: item[2]
            ]
        }

        wgs_hap2 = HAPLOTYPE_READS_WGS.out.hap2.map { item ->
            [
                name: item[0],
                r1: item[1],
                r2: item[2]
            ]
        }


        // -------------------------------------------------------------------
        // HAP1
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi_hap1 = ALIGN_HIFI(
                ready,
                hap1,
                hifi_hap1
            )
        }

        if (params.ont) {

            align_ont_hap1 = ALIGN_ONT(
                ready,
                hap1,
                ont_hap1
            )
        }

        align_wgs_hap1 = ALIGN_WGS(
            ready,
            hap1,
            wgs_hap1.r1,
            wgs_hap1.r2
        )


        // -------------------------------------------------------------------
        // HAP2
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi_hap2 = ALIGN_HIFI(
                ready,
                hap2,
                hifi_hap2
            )
        }

        if (params.ont) {

            align_ont_hap2 = ALIGN_ONT(
                ready,
                hap2,
                ont_hap2
            )
        }

        align_wgs_hap2 = ALIGN_WGS(
            ready,
            hap2,
            wgs_hap2.r1,
            wgs_hap2.r2
        )


        // -------------------------------------------------------------------
        // DeepVariant
        //
        // Отдельный DeepVariant для каждого haplotype.
        // -------------------------------------------------------------------

        if (params.hifi) {

            deepvariant_hap1 = DEEPVARIANT(
                ready,
                hap1,
                hap1_fai,
                align_wgs_hap1.bam,
                align_wgs_hap1.bai,
                align_hifi_hap1.bam
            )

            deepvariant_hap2 = DEEPVARIANT(
                ready,
                hap2,
                hap2_fai,
                align_wgs_hap2.bam,
                align_wgs_hap2.bai,
                align_hifi_hap2.bam
            )

        } else {

            deepvariant_hap1 = DEEPVARIANT(
                ready,
                hap1,
                hap1_fai,
                align_wgs_hap1.bam,
                align_wgs_hap1.bai,
                "none"
            )

            deepvariant_hap2 = DEEPVARIANT(
                ready,
                hap2,
                hap2_fai,
                align_wgs_hap2.bam,
                align_wgs_hap2.bai,
                "none"
            )
        }


        // -------------------------------------------------------------------
        // PEPPER
        //
        // Только если есть ONT.
        // -------------------------------------------------------------------

        if (params.ont) {

            pepper_hap1 = PEPPER(
                ready,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                hap1,
                hap1_fai
            )

            pepper_hap2 = PEPPER(
                ready,
                align_ont_hap2.bam,
                align_ont_hap2.bai,
                hap2,
                hap2_fai
            )
        }


        // -------------------------------------------------------------------
        // ORIGINAL T2T
        // -------------------------------------------------------------------

        if (params.ont) {

            original_t2t_hap1 = ORIGINAL_T2T(
                ready,
                deepvariant_hap1.vcf,
                pepper_hap1.vcf,
                hap1,
                wgs_hap1.r1,
                wgs_hap1.r2,
                params.hifi ? hifi_hap1 : "none"
            )

            original_t2t_hap2 = ORIGINAL_T2T(
                ready,
                deepvariant_hap2.vcf,
                pepper_hap2.vcf,
                hap2,
                wgs_hap2.r1,
                wgs_hap2.r2,
                params.hifi ? hifi_hap2 : "none"
            )

        } else {

            original_t2t_hap1 = ORIGINAL_T2T(
                ready,
                deepvariant_hap1.vcf,
                "none",
                hap1,
                wgs_hap1.r1,
                wgs_hap1.r2,
                params.hifi ? hifi_hap1 : "none"
            )

            original_t2t_hap2 = ORIGINAL_T2T(
                ready,
                deepvariant_hap2.vcf,
                "none",
                hap2,
                wgs_hap2.r1,
                wgs_hap2.r2,
                params.hifi ? hifi_hap2 : "none"
            )
        }


        // -------------------------------------------------------------------
        // AUTO POLISH
        // -------------------------------------------------------------------

        if (params.hifi) {

            auto_polish_hap1 = AUTO_POLISH(
                ready,
                hap1,
                hifi_hap1,
                original_t2t_hap1.readmers_meryl,
                "hifi"
            )

            auto_polish_hap2 = AUTO_POLISH(
                ready,
                hap2,
                hifi_hap2,
                original_t2t_hap2.readmers_meryl,
                "hifi"
            )

        } else if (params.ont) {

            auto_polish_hap1 = AUTO_POLISH(
                ready,
                hap1,
                ont_hap1,
                original_t2t_hap1.readmers_meryl,
                "ont"
            )

            auto_polish_hap2 = AUTO_POLISH(
                ready,
                hap2,
                ont_hap2,
                original_t2t_hap2.readmers_meryl,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // NextPolish2
        //
        // Сохраняем исходное разделение:
        // HiFi есть  -> ALIGN_HIFI
        // HiFi нет   -> ALIGN_ONT
        // -------------------------------------------------------------------

        if (params.hifi) {

            np2_hap1 = NP2(
                ready,
                hap1,
                align_hifi_hap1.bam,
                wgs_hap1.r1,
                wgs_hap1.r2
            )

            np2_hap2 = NP2(
                ready,
                hap2,
                align_hifi_hap2.bam,
                wgs_hap2.r1,
                wgs_hap2.r2
            )

        } else {

            np2_hap1 = NP2(
                ready,
                hap1,
                align_ont_hap1.bam,
                wgs_hap1.r1,
                wgs_hap1.r2
            )

            np2_hap2 = NP2(
                ready,
                hap2,
                align_ont_hap2.bam,
                wgs_hap2.r1,
                wgs_hap2.r2
            )
        }


        // -------------------------------------------------------------------
        // Sniffles
        // -------------------------------------------------------------------

        if (params.hifi) {

            sniffles_hap1 = SNIFFLES_HIFI(
                ready,
                hap1,
                align_hifi_hap1.bam,
                align_hifi_hap1.bai,
                "hifi"
            )

            sniffles_hap2 = SNIFFLES_HIFI(
                ready,
                hap2,
                align_hifi_hap2.bam,
                align_hifi_hap2.bai,
                "hifi"
            )

        } else {

            sniffles_hap1 = SNIFFLES_ONT(
                ready,
                hap1,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                "ont"
            )

            sniffles_hap2 = SNIFFLES_ONT(
                ready,
                hap2,
                align_ont_hap2.bam,
                align_ont_hap2.bai,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // CuteSV
        // -------------------------------------------------------------------

        if (params.hifi) {

            cutesv_hap1 = CUTESV_HIFI(
                ready,
                hap1,
                align_hifi_hap1.bam,
                align_hifi_hap1.bai,
                "hifi"
            )

            cutesv_hap2 = CUTESV_HIFI(
                ready,
                hap2,
                align_hifi_hap2.bam,
                align_hifi_hap2.bai,
                "hifi"
            )

        } else {

            cutesv_hap1 = CUTESV_ONT(
                ready,
                hap1,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                "ont"
            )

            cutesv_hap2 = CUTESV_ONT(
                ready,
                hap2,
                align_ont_hap2.bam,
                align_ont_hap2.bai,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // Flagger
        //
        // HiFi -> HiFi BAM
        // no HiFi -> ONT BAM
        // -------------------------------------------------------------------

        if (params.hifi) {

            flagger_hap1 = FLAGGER(
                ready,
                hap1,
                align_hifi_hap1.bam
            )

            flagger_hap2 = FLAGGER(
                ready,
                hap2,
                align_hifi_hap2.bam
            )

        } else {

            flagger_hap1 = FLAGGER(
                ready,
                hap1,
                align_ont_hap1.bam
            )

            flagger_hap2 = FLAGGER(
                ready,
                hap2,
                align_ont_hap2.bam
            )
        }


        // -------------------------------------------------------------------
        // Merqury
        // -------------------------------------------------------------------

        merqury_hap1 = MERQURY(
            ready,
            hap1,
            wgs_hap1.r1,
            wgs_hap1.r2
        )

        merqury_hap2 = MERQURY(
            ready,
            hap2,
            wgs_hap2.r1,
            wgs_hap2.r2
        )


        // -------------------------------------------------------------------
        // Medaka
        //
        // Только при наличии ONT.
        // -------------------------------------------------------------------

        if (params.ont) {

            medaka_hap1 = MEDAKA(
                ready,
                hap1,
                ont_hap1
            )

            medaka_hap2 = MEDAKA(
                ready,
                hap2,
                ont_hap2
            )
        }


    } else {


        // ===================================================================
        // HAPLOID — оригинальная ветка
        // ===================================================================

        draft = Channel.value(
            file(params.draft, checkIfExists: true)
        )

        draft_fai = Channel.value(
            file("${params.draft}.fai", checkIfExists: true)
        )


        // -------------------------------------------------------------------
        // Alignments
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi = ALIGN_HIFI(
                ready,
                draft,
                params.hifi
            )
        }

        if (params.ont) {

            align_ont = ALIGN_ONT(
                ready,
                draft,
                params.ont
            )
        }

        align_wgs = ALIGN_WGS(
            ready,
            draft,
            wgs_r1,
            wgs_r2
        )


        // -------------------------------------------------------------------
        // DeepVariant
        // -------------------------------------------------------------------

        if (params.hifi) {

            deepvariant = DEEPVARIANT(
                ready,
                draft,
                draft_fai,
                align_wgs.bam,
                align_wgs.bai,
                align_hifi.bam
            )

        } else {

            deepvariant = DEEPVARIANT(
                ready,
                draft,
                draft_fai,
                align_wgs.bam,
                align_wgs.bai,
                "none"
            )
        }


        // -------------------------------------------------------------------
        // Pepper
        // -------------------------------------------------------------------

        if (params.ont) {

            pepper = PEPPER(
                ready,
                align_ont.bam,
                align_ont.bai,
                draft,
                draft_fai
            )
        }


        // -------------------------------------------------------------------
        // ORIGINAL T2T
        // -------------------------------------------------------------------

        if (params.ont) {

            original_t2t = ORIGINAL_T2T(
                ready,
                deepvariant.vcf,
                pepper.vcf,
                draft,
                wgs_r1,
                wgs_r2,
                params.hifi ? params.hifi : "none"
            )

        } else {

            original_t2t = ORIGINAL_T2T(
                ready,
                deepvariant.vcf,
                "none",
                draft,
                wgs_r1,
                wgs_r2,
                params.hifi ? params.hifi : "none"
            )
        }


        // -------------------------------------------------------------------
        // Auto polish
        // -------------------------------------------------------------------

        if (params.hifi) {

            auto_polish = AUTO_POLISH(
                ready,
                draft,
                params.hifi,
                original_t2t.readmers_meryl,
                "hifi"
            )

        } else if (params.ont) {

            auto_polish = AUTO_POLISH(
                ready,
                draft,
                params.ont,
                original_t2t.readmers_meryl,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // NextPolish2
        // -------------------------------------------------------------------

        if (params.hifi) {

            np2 = NP2(
                ready,
                draft,
                align_hifi.bam,
                wgs_r1,
                wgs_r2
            )

        } else {

            np2 = NP2(
                ready,
                draft,
                align_ont.bam,
                wgs_r1,
                wgs_r2
            )
        }


        // -------------------------------------------------------------------
        // Sniffles
        // -------------------------------------------------------------------

        if (params.hifi) {

            sniffles_hifi = SNIFFLES_HIFI(
                ready,
                draft,
                align_hifi.bam,
                align_hifi.bai,
                "hifi"
            )

        } else {

            sniffles_ont = SNIFFLES_ONT(
                ready,
                draft,
                align_ont.bam,
                align_ont.bai,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // CuteSV
        // -------------------------------------------------------------------

        if (params.hifi) {

            cutesv_hifi = CUTESV_HIFI(
                ready,
                draft,
                align_hifi.bam,
                align_hifi.bai,
                "hifi"
            )

        } else {

            cutesv_ont = CUTESV_ONT(
                ready,
                draft,
                align_ont.bam,
                align_ont.bai,
                "ont"
            )
        }


        // -------------------------------------------------------------------
        // Flagger
        // -------------------------------------------------------------------

        if (params.hifi) {

            flagger = FLAGGER(
                ready,
                draft,
                align_hifi.bam
            )

        } else {

            flagger = FLAGGER(
                ready,
                draft,
                align_ont.bam
            )
        }


        // -------------------------------------------------------------------
        // Merqury
        // -------------------------------------------------------------------

        merqury = MERQURY(
            ready,
            draft,
            wgs_r1,
            wgs_r2
        )


        // -------------------------------------------------------------------
        // Medaka
        // -------------------------------------------------------------------

        if (params.ont) {

            medaka = MEDAKA(
                ready,
                draft,
                params.ont
            )
        }
    }
}