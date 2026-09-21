nextflow.enable.dsl = 2

include { CHECK        } from './modules/check_containers'
include { ALIGN_HIFI } from './modules/align'
include { ALIGN_ONT }  from './modules/align'
include { ALIGN_ILLUMINA } from './modules/align'
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

workflow {
    ready       = CHECK()
    draft       = Channel.fromPath(params.draft)
    draft_fai = Channel.fromPath(params.draft+".fai")
    ont         = Channel.fromPath(params.ont)

    illumina_r1 = Channel.fromPath(params.illumina_r1)
    illumina_r2 = Channel.fromPath(params.illumina_r2)

    align_illumina = ALIGN_ILLUMINA( 
        ready, 
        draft, 
        illumina_r1, 
        illumina_r2 
    )

    align_ont = ALIGN_ONT( 
        ready, 
        draft, 
        ont 
    )

    if (params.hifi) {
        
        hifi = Channel.fromPath( 
            params.hifi, 
            checkIfExists: true 
        )

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
            illumina_r1,
            illumina_r2,
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
            illumina_r1,
            illumina_r2
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
            align_hifi.bai,
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
            align_hifi.bai
        )

        merqury = MERQURY(
            ready,
            draft,
            illumina_r1,
            illumina_r2
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
            illumina_r1,
            illumina_r2
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
            illumina_r1,
            illumina_r2
        )

        t2t = ORIGINAL_T2T(
            ready,
            deepvariant.vcf,
            pepper.vcf,
            draft,
            illumina_r1,
            illumina_r2,
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

