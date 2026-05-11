nextflow.enable.dsl = 2

include { CHECK        } from './modules/check_containers'
include { ALIGN_ALL        } from './modules/align'
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

workflow {
    ready       = CHECK()
    draft       = Channel.fromPath(params.draft)
    draft_fai = Channel.fromPath(params.draft+".fai")
    hifi        = Channel.fromPath(params.hifi)
    ont         = Channel.fromPath(params.ont)

    illumina_r1 = Channel.fromPath(params.illumina_r1)
    illumina_r2 = Channel.fromPath(params.illumina_r2)

    align = ALIGN_ALL(
        ready,
        draft,
        hifi,
        ont,
        illumina_r1,
        illumina_r2
    )

    deepvariant = DEEPVARIANT(
        ready,
        draft,
        draft_fai,
        align.hifi_bam,
        align.illumina_bam
    )

    pepper = PEPPER(
        ready,
        align.ont_bam,
        align.ont_bai,
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
        t2t.readmers_meryl
    )

    np2 = NP2(
        ready,
        draft,
        align.hifi_bam,
        illumina_r1,
        illumina_r2
    )

    sniffles_ont = SNIFFLES_ONT(
        ready,
        draft,
        align.ont_bam,
        align.ont_bai,
        "ont"
    )

    sniffles_hifi = SNIFFLES_HIFI(
        ready,
        draft,
        align.hifi_bam,
        align.hifi_bai,
        "hifi"
    )

    cutesv_ont = CUTESV_ONT(
        ready,
        draft,
        align.ont_bam,
        align.ont_bai,
        "ont"
    )

    cutesv_hifi = CUTESV_HIFI(
        ready,
        draft,
        align.hifi_bam,
        align.hifi_bai,
        "hifi"
    )

    flagger = FLAGGER(
        ready,
        draft,
        align.hifi_bam
    )

    merqury = MERQURY(
        ready,
        draft,
        illumina_r1,
        illumina_r2
    )

}

