nextflow.enable.dsl = 2

include { ALIGN } from './modules/align'
include { PARLIAMENT2 } from './modules/parliament2'
include { SNIFFLES } from './modules/sniffles'
include { DEEPVARIANT } from './modules/deepvariant'
include { PEPPER } from './modules/pepper'
include { JASMINE } from './modules/jasmine'
include { MERFIN } from './modules/merfin'
include { IRIS } from './modules/iris'

include { REPEATMASKER } from './modules/repeatmasker'
include { FLAGGER } from './modules/flagger'

include { CATBOOST } from './modules/catboost'


workflow {

    assembly_ch = Channel.fromPath(params.assembly)
    reads_ch = Channel
        .fromPath(params.reads_table)
        .splitCsv(header:true, sep:/\s+/)
        .map { row -> tuple(file(row.reads_path), row.flag) }

    /*
    --------------------
    ALIGN
    --------------------
    */

    align_out = ALIGN(assembly_ch, reads_ch)

    illumina_bam = align_out.filter{f,b -> f=="illumina"}.map{f,b->b}
    hifi_bam     = align_out.filter{f,b -> f=="hifi"}.map{f,b->b}
    ont_bam      = align_out.filter{f,b -> f=="ont"}.map{f,b->b}
    clr_bam      = align_out.filter{f,b -> f=="clr"}.map{f,b->b}


    /*
    --------------------
    SV CALLERS
    --------------------
    */

    parliament_vcf = PARLIAMENT2(illumina_bam)

    sniffles_vcf = SNIFFLES(
        Channel.merge(hifi_bam, ont_bam, clr_bam)
    )

    sv_vcfs = Channel.merge(
        parliament_vcf,
        sniffles_vcf
    )


    /*
    --------------------
    SNV CALLERS
    --------------------
    */

    deepvariant_vcf = DEEPVARIANT(
        Channel.merge(illumina_bam, hifi_bam)
    )

    pepper_vcf = PEPPER(ont_bam)

    snv_vcfs = Channel.merge(
        deepvariant_vcf,
        pepper_vcf
    )


    /*
    --------------------
    MERGE VARIANTS
    --------------------
    */

    sv_merged = JASMINE(sv_vcfs.collect())
    snv_merged = JASMINE(snv_vcfs.collect())


    /*
    --------------------
    MERFIN
    --------------------
    */

    merfin_out = MERFIN(
        assembly_ch,
        sv_merged
    )


    /*
    --------------------
    IRIS (SV polishing)
    --------------------
    */

    IRIS(
        assembly_ch,
        hifi_bam.first(),
        sv_merged,
        parliament_vcf
    )


    /*
    --------------------
    GENOME ANNOTATION
    --------------------
    */

    repeat_gff = REPEATMASKER(assembly_ch)

    flagger_out = FLAGGER(
        assembly_ch,
        hifi_bam.first()
    )


    /*
    --------------------
    CATBOOST
    --------------------
    */

    CATBOOST(
        Channel.value("SV"),
        sv_merged,
        repeat_gff,
        flagger_out,
        merfin_out
    )

    CATBOOST(
        Channel.value("SNV"),
        snv_merged,
        repeat_gff,
        flagger_out,
        merfin_out
    )

}