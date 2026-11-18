nextflow.enable.dsl = 2

include { CHECK        } from './modules/check_containers'

// ---------------------------------------------------------------------------
// Haploid / общие имена (каждый процесс вызывается максимум один раз)
// ---------------------------------------------------------------------------
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
include { ML_PREDICT } from './modules/ml_predict'
include { ML_TRAIN } from './modules/ml_train'
include { TRUTH_VCF } from './modules/truth_vcf'
include { EMPTY_VCF } from './modules/empty_vcf'

// ---------------------------------------------------------------------------
// Diploid: отдельный алиас на каждый haplotype (H1 / H2)
// ---------------------------------------------------------------------------
include { ALIGN_HIFI as ALIGN_HIFI_H1 } from './modules/align'
include { ALIGN_HIFI as ALIGN_HIFI_H2 } from './modules/align'
include { ALIGN_ONT  as ALIGN_ONT_H1  } from './modules/align'
include { ALIGN_ONT  as ALIGN_ONT_H2  } from './modules/align'
include { ALIGN_WGS  as ALIGN_WGS_H1  } from './modules/align'
include { ALIGN_WGS  as ALIGN_WGS_H2  } from './modules/align'

include { DEEPVARIANT as DEEPVARIANT_H1 } from './modules/deepvariant'
include { DEEPVARIANT as DEEPVARIANT_H2 } from './modules/deepvariant'

include { PEPPER as PEPPER_H1 } from './modules/pepper'
include { PEPPER as PEPPER_H2 } from './modules/pepper'

include { ORIGINAL_T2T as ORIGINAL_T2T_H1 } from './modules/original_t2t'
include { ORIGINAL_T2T as ORIGINAL_T2T_H2 } from './modules/original_t2t'

include { AUTO_POLISH as AUTO_POLISH_H1 } from './modules/automated_polishing'
include { AUTO_POLISH as AUTO_POLISH_H2 } from './modules/automated_polishing'

include { NP2 as NP2_H1 } from './modules/nextpolish2'
include { NP2 as NP2_H2 } from './modules/nextpolish2'

include { SNIFFLES as SNIFFLES_HIFI_H1 } from './modules/sniffles'
include { SNIFFLES as SNIFFLES_HIFI_H2 } from './modules/sniffles'
include { SNIFFLES as SNIFFLES_ONT_H1  } from './modules/sniffles'
include { SNIFFLES as SNIFFLES_ONT_H2  } from './modules/sniffles'

include { CUTESV as CUTESV_HIFI_H1 } from './modules/cutesv'
include { CUTESV as CUTESV_HIFI_H2 } from './modules/cutesv'
include { CUTESV as CUTESV_ONT_H1  } from './modules/cutesv'
include { CUTESV as CUTESV_ONT_H2  } from './modules/cutesv'

include { FLAGGER as FLAGGER_H1 } from './modules/flagger'
include { FLAGGER as FLAGGER_H2 } from './modules/flagger'

include { MERQURY as MERQURY_H1 } from './modules/merqury'
include { MERQURY as MERQURY_H2 } from './modules/merqury'

include { MEDAKA as MEDAKA_H1 } from './modules/medaka'
include { MEDAKA as MEDAKA_H2 } from './modules/medaka'
include { ML_PREDICT as ML_PREDICT_H1 } from './modules/ml_predict'
include { ML_PREDICT as ML_PREDICT_H2 } from './modules/ml_predict'
include { ML_TRAIN as ML_TRAIN_H1 } from './modules/ml_train'
include { ML_TRAIN as ML_TRAIN_H2 } from './modules/ml_train'
include { TRUTH_VCF as TRUTH_VCF_H1 } from './modules/truth_vcf'
include { TRUTH_VCF as TRUTH_VCF_H2 } from './modules/truth_vcf'
include { EMPTY_VCF as EMPTY_VCF_H1 } from './modules/empty_vcf'
include { EMPTY_VCF as EMPTY_VCF_H2 } from './modules/empty_vcf'

// ---------------------------------------------------------------------------
// Уже раздельные процессы
// ---------------------------------------------------------------------------
include { MERGE_READS as MERGE_WGS_R1 } from './modules/merge_reads'
include { MERGE_READS as MERGE_WGS_R2 } from './modules/merge_reads'
include { HAPLOTYPE_READS_HIFI } from './modules/haplotype_reads'
include { HAPLOTYPE_READS_ONT  } from './modules/haplotype_reads'
include { HAPLOTYPE_READS_WGS  } from './modules/haplotype_reads'


// Принимает файл, glob или ДИРЕКТОРИЮ с ридами и возвращает канал файлов.
def reads_channel(p) {
    def inputs = p instanceof Collection ? p.flatten() : [p]
    def channels = inputs.collect { rawInput ->
        def inputPath = rawInput.toString()

        // Avoid calling isDirectory() on a glob, and accept either one path
        // or a list of files/directories/globs from the params file.
        if (!(inputPath =~ /[*?{}]/)) {
            def candidate = file(inputPath)
            if (candidate.isDirectory()) {
                return Channel.fromPath(
                    "${candidate}/**/*.{fastq,fq,fastq.gz,fq.gz}",
                    checkIfExists: true
                )
            }
        }

        Channel.fromPath(inputPath, checkIfExists: true)
    }
    mix_channels(channels)
}

def as_bool(value) {
    value instanceof Boolean ? value : value?.toString()?.equalsIgnoreCase('true')
}

def mix_channels(channels) {
    def active = channels.findAll { it != null }
    if (active.isEmpty()) return Channel.empty()
    active.drop(1).inject(active[0]) { merged, next -> merged.mix(next) }
}

workflow {

    run_deepvariant = as_bool(params.run_deepvariant)
    run_pepper = as_bool(params.run_pepper)
    run_original_t2t = as_bool(params.run_original_t2t)
    run_auto_polish = as_bool(params.run_auto_polish)
    run_nextpolish2 = as_bool(params.run_nextpolish2)
    run_sniffles = as_bool(params.run_sniffles)
    run_cutesv = as_bool(params.run_cutesv)
    run_flagger = as_bool(params.run_flagger)
    run_merqury = as_bool(params.run_merqury)
    run_medaka = as_bool(params.run_medaka)

    def read_mode = params.hifi && params.ont ? 'WGS+HiFi+ONT' :
                    params.hifi ? 'WGS+HiFi' : 'WGS+ONT'
    log.info "Read mode: ${read_mode}"

    if (!params.hifi && !params.ont) {
        error "ERROR: provide at least one long-read input: --hifi and/or --ont"
    }
    workflow_mode = params.workflow_mode?.toString()?.toLowerCase()
    if (!(workflow_mode in ['predict', 'train'])) {
        error "ERROR: --workflow_mode must be predict or train"
    }
    if (!(params.ploidy in ['haploid', 'diploid'])) {
        error "ERROR: --ploidy must be haploid or diploid"
    }
    def ml_threshold = params.ml_threshold.toString().toDouble()
    if (workflow_mode == 'predict' && (ml_threshold < 0 || ml_threshold > 1)) {
        error "ERROR: --ml_threshold must be between 0 and 1"
    }
    if (workflow_mode == 'predict') {
        model_paths = params.ploidy == 'diploid' ?
            [params.ml_model_hap1 ?: params.ml_model, params.ml_model_hap2 ?: params.ml_model] :
            [params.ml_model]
        if (model_paths.any { !it }) {
            error "ERROR: prediction requires --ml_model (or both haplotype-specific model paths)"
        }
        if (model_paths.any { !it.toString().toLowerCase().endsWith('.cbm') }) {
            error "ERROR: model paths must point to CatBoost .cbm files"
        }
    } else if (params.ploidy == 'haploid') {
        if (!params.truth_fasta) error "ERROR: train mode requires --truth_fasta"
    } else if (!params.truth_fasta_hap1 || !params.truth_fasta_hap2) {
        error "ERROR: diploid train mode requires --truth_fasta_hap1 and --truth_fasta_hap2"
    }
    if (run_original_t2t && !run_deepvariant) {
        error "ERROR: --run_original_t2t requires --run_deepvariant true"
    }
    if (!(run_deepvariant || (run_pepper && params.ont) || run_original_t2t ||
          (run_nextpolish2 && params.hifi) || run_sniffles || run_cutesv)) {
        error "ERROR: enable at least one VCF-producing caller/polisher for ML input"
    }

    containers_to_check = [params.mapping_container, params.ml_container]
    if (run_deepvariant) containers_to_check << params.deepvariant_container
    if (run_pepper && params.ont) containers_to_check << params.pepper_container
    if (run_original_t2t || run_auto_polish ||
        (run_nextpolish2 && params.hifi) || (run_medaka && params.ont) || run_merqury) {
        containers_to_check << params.polishing_container
    }
    if (run_sniffles || run_cutesv) containers_to_check << params.sv_container
    if (run_flagger) containers_to_check << params.flagger_container
    if (workflow_mode == 'train') containers_to_check << params.truth_container

    ready = CHECK(containers_to_check.unique())
    ml_script = Channel.value(file("${baseDir}/scripts/vcf_catboost.py", checkIfExists: true))
    truth_script = Channel.value(file("${baseDir}/../my_scripts/truth_paftools_syri.sh", checkIfExists: true))
    if (workflow_mode == 'predict') {
        if (params.ploidy == 'diploid') {
            model_path_hap1 = params.ml_model_hap1 ?: params.ml_model
            model_path_hap2 = params.ml_model_hap2 ?: params.ml_model
            ml_model_hap1 = Channel.value(file(model_path_hap1, checkIfExists: true))
            ml_model_hap2 = Channel.value(file(model_path_hap2, checkIfExists: true))
        } else {
            ml_model = Channel.value(file(params.ml_model, checkIfExists: true))
        }
    }


    // -----------------------------------------------------------------------
    // WGS reads
    // -----------------------------------------------------------------------

    wgs_path = params.wgs ?: params.illumina ?: params.mgi
    if (!wgs_path) {
        error "ERROR: provide paired WGS reads with --wgs (or legacy --illumina/--mgi)"
    }

    wgs_files = reads_channel(wgs_path).collect()
    wgs_pairs = wgs_files.map { files ->
        def sorted = files.sort { it.name }
        def fastqFiles = sorted.findAll { it.name ==~ /(?i).+\.(?:fastq|fq)(?:\.gz)?$/ }
        def r1 = sorted.findAll { it.name ==~ /(?i).+_(?:R1|1)(?:_[0-9]+)?\.(?:fastq|fq)(?:\.gz)?$/ }
        def r2 = sorted.findAll { it.name ==~ /(?i).+_(?:R2|2)(?:_[0-9]+)?\.(?:fastq|fq)(?:\.gz)?$/ }
        def keyFor = { f, readNumber ->
            def stem = f.name.replaceFirst(/(?i)\.(?:fastq|fq)(?:\.gz)?$/, '')
            def markerPattern = readNumber == 1 ? '(?i)_(?:R1|1)(_[0-9]+)?$' : '(?i)_(?:R2|2)(_[0-9]+)?$'
            def markerPatternCompiled = java.util.regex.Pattern.compile(markerPattern)
            def marker = markerPatternCompiled.matcher(stem)
            if (!marker.find()) return null
            def lane = marker.group(1) ?: ''
            def sample = stem.substring(0, marker.start())
            lane ? "${sample}${lane}" : sample
        }
        def r1ByKey = [:]
        def r2ByKey = [:]
        r1.each { f ->
            def key = keyFor(f, 1)
            if (r1ByKey.containsKey(key)) throw new IllegalArgumentException("Duplicate WGS R1 sample/lane: ${key}")
            r1ByKey[key] = f
        }
        r2.each { f ->
            def key = keyFor(f, 2)
            if (r2ByKey.containsKey(key)) throw new IllegalArgumentException("Duplicate WGS R2 sample/lane: ${key}")
            r2ByKey[key] = f
        }
        if (r1ByKey.isEmpty() || r1ByKey.keySet() != r2ByKey.keySet() ||
            fastqFiles.size() != r1.size() + r2.size()) {
            throw new IllegalArgumentException("WGS input must contain matching R1/R2 FASTQ files with _R1/_R2 or _1/_2 suffixes")
        }
        def keys = r1ByKey.keySet().sort()
        tuple(keys.collect { r1ByKey[it] }, keys.collect { r2ByKey[it] })
    }

    wgs_r1_raw = wgs_pairs.map { it[0] }
    wgs_r2_raw = wgs_pairs.map { it[1] }


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

        if (workflow_mode == 'train') {
            truth_query_hap1 = Channel.value(file(params.truth_fasta_hap1, checkIfExists: true))
            truth_query_hap2 = Channel.value(file(params.truth_fasta_hap2, checkIfExists: true))
            truth_vcf_hap1 = TRUTH_VCF_H1(ready, truth_script, hap1, truth_query_hap1, "hap1").vcf
            truth_vcf_hap2 = TRUTH_VCF_H2(ready, truth_script, hap2, truth_query_hap2, "hap2").vcf
        }

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

            hifi_reads = reads_channel(params.hifi).collect()

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

            ont_reads = reads_channel(params.ont).collect()

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

        // ИСПРАВЛЕНО: раньше было wgs_hap1.r1 на канале (так нельзя) —
        // теперь отдельные каналы R1 / R2 для каждого haplotype.

        wgs_hap1_r1 = HAPLOTYPE_READS_WGS.out.hap1.map { it[1] }
        wgs_hap1_r2 = HAPLOTYPE_READS_WGS.out.hap1.map { it[2] }

        wgs_hap2_r1 = HAPLOTYPE_READS_WGS.out.hap2.map { it[1] }
        wgs_hap2_r2 = HAPLOTYPE_READS_WGS.out.hap2.map { it[2] }


        // -------------------------------------------------------------------
        // HAP1
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi_hap1 = ALIGN_HIFI_H1(
                ready,
                hap1,
                hifi_hap1
            )
        }

        if (params.ont) {

            align_ont_hap1 = ALIGN_ONT_H1(
                ready,
                hap1,
                ont_hap1
            )
        }

        align_wgs_hap1 = ALIGN_WGS_H1(
            ready,
            hap1,
            wgs_hap1_r1,
            wgs_hap1_r2
        )


        // -------------------------------------------------------------------
        // HAP2
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi_hap2 = ALIGN_HIFI_H2(
                ready,
                hap2,
                hifi_hap2
            )
        }

        if (params.ont) {

            align_ont_hap2 = ALIGN_ONT_H2(
                ready,
                hap2,
                ont_hap2
            )
        }

        align_wgs_hap2 = ALIGN_WGS_H2(
            ready,
            hap2,
            wgs_hap2_r1,
            wgs_hap2_r2
        )


        // -------------------------------------------------------------------
        // DeepVariant
        //
        // Отдельный DeepVariant для каждого haplotype.
        // -------------------------------------------------------------------

        if (run_deepvariant && params.hifi) {

            deepvariant_hap1 = DEEPVARIANT_H1(
                ready,
                hap1,
                hap1_fai,
                align_wgs_hap1.bam,
                align_wgs_hap1.bai,
                align_hifi_hap1.bam,
                true
            )

            deepvariant_hap2 = DEEPVARIANT_H2(
                ready,
                hap2,
                hap2_fai,
                align_wgs_hap2.bam,
                align_wgs_hap2.bai,
                align_hifi_hap2.bam,
                true
            )

        } else if (run_deepvariant) {

            deepvariant_hap1 = DEEPVARIANT_H1(
                ready,
                hap1,
                hap1_fai,
                align_wgs_hap1.bam,
                align_wgs_hap1.bai,
                align_wgs_hap1.bam,
                false
            )

            deepvariant_hap2 = DEEPVARIANT_H2(
                ready,
                hap2,
                hap2_fai,
                align_wgs_hap2.bam,
                align_wgs_hap2.bai,
                align_wgs_hap2.bam,
                false
            )
        }


        // -------------------------------------------------------------------
        // PEPPER
        //
        // Только если есть ONT.
        // -------------------------------------------------------------------

        if (run_pepper && params.ont) {

            pepper_hap1 = PEPPER_H1(
                ready,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                hap1,
                hap1_fai
            )

            pepper_hap2 = PEPPER_H2(
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

        if (run_original_t2t) {
            if (run_pepper && params.ont) {
                original_pepper_hap1 = pepper_hap1.vcf
                original_pepper_hap2 = pepper_hap2.vcf
            } else {
                empty_vcf_hap1 = EMPTY_VCF_H1(deepvariant_hap1.vcf)
                empty_vcf_hap2 = EMPTY_VCF_H2(deepvariant_hap2.vcf)
                original_pepper_hap1 = empty_vcf_hap1.vcf
                original_pepper_hap2 = empty_vcf_hap2.vcf
            }

            original_t2t_hap1 = ORIGINAL_T2T_H1(
                ready,
                deepvariant_hap1.vcf,
                original_pepper_hap1,
                hap1,
                wgs_hap1_r1,
                wgs_hap1_r2,
                params.hifi ? hifi_hap1 : wgs_hap1_r1,
                params.hifi ? true : false
            )

            original_t2t_hap2 = ORIGINAL_T2T_H2(
                ready,
                deepvariant_hap2.vcf,
                original_pepper_hap2,
                hap2,
                wgs_hap2_r1,
                wgs_hap2_r2,
                params.hifi ? hifi_hap2 : wgs_hap2_r1,
                params.hifi ? true : false
            )
        }


        // -------------------------------------------------------------------
        // AUTO POLISH
        // -------------------------------------------------------------------

        if (run_auto_polish) {
            if (params.hifi) {
                auto_polish_hap1 = AUTO_POLISH_H1(
                    ready, hap1, hifi_hap1, wgs_hap1_r1, wgs_hap1_r2, "hifi"
                )
                auto_polish_hap2 = AUTO_POLISH_H2(
                    ready, hap2, hifi_hap2, wgs_hap2_r1, wgs_hap2_r2, "hifi"
                )
            } else {
                auto_polish_hap1 = AUTO_POLISH_H1(
                    ready, hap1, ont_hap1, wgs_hap1_r1, wgs_hap1_r2, "ont"
                )
                auto_polish_hap2 = AUTO_POLISH_H2(
                    ready, hap2, ont_hap2, wgs_hap2_r1, wgs_hap2_r2, "ont"
                )
            }
        }


        // -------------------------------------------------------------------
        // NextPolish2
        //
        // Сохраняем исходное разделение:
        // HiFi есть  -> ALIGN_HIFI
        // HiFi нет   -> ALIGN_ONT
        // -------------------------------------------------------------------

        if (run_nextpolish2 && params.hifi) {

            np2_hap1 = NP2_H1(
                ready,
                hap1,
                align_hifi_hap1.bam,
                wgs_hap1_r1,
                wgs_hap1_r2
            )

            np2_hap2 = NP2_H2(
                ready,
                hap2,
                align_hifi_hap2.bam,
                wgs_hap2_r1,
                wgs_hap2_r2
            )

        }


        // -------------------------------------------------------------------
        // Sniffles
        // -------------------------------------------------------------------

        if (run_sniffles && params.hifi) {

            sniffles_hifi_hap1 = SNIFFLES_HIFI_H1(
                ready,
                hap1,
                align_hifi_hap1.bam,
                align_hifi_hap1.bai,
                "hifi"
            )

            sniffles_hifi_hap2 = SNIFFLES_HIFI_H2(
                ready,
                hap2,
                align_hifi_hap2.bam,
                align_hifi_hap2.bai,
                "hifi"
            )

        }

        if (run_sniffles && params.ont) {

            sniffles_ont_hap1 = SNIFFLES_ONT_H1(
                ready,
                hap1,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                "ont"
            )

            sniffles_ont_hap2 = SNIFFLES_ONT_H2(
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

        if (run_cutesv && params.hifi) {

            cutesv_hifi_hap1 = CUTESV_HIFI_H1(
                ready,
                hap1,
                align_hifi_hap1.bam,
                align_hifi_hap1.bai,
                "hifi"
            )

            cutesv_hifi_hap2 = CUTESV_HIFI_H2(
                ready,
                hap2,
                align_hifi_hap2.bam,
                align_hifi_hap2.bai,
                "hifi"
            )

        }

        if (run_cutesv && params.ont) {

            cutesv_ont_hap1 = CUTESV_ONT_H1(
                ready,
                hap1,
                align_ont_hap1.bam,
                align_ont_hap1.bai,
                "ont"
            )

            cutesv_ont_hap2 = CUTESV_ONT_H2(
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

        if (run_flagger && params.hifi) {

            flagger_hap1 = FLAGGER_H1(
                ready,
                hap1,
                align_hifi_hap1.bam
            )

            flagger_hap2 = FLAGGER_H2(
                ready,
                hap2,
                align_hifi_hap2.bam
            )

        } else if (run_flagger) {

            flagger_hap1 = FLAGGER_H1(
                ready,
                hap1,
                align_ont_hap1.bam
            )

            flagger_hap2 = FLAGGER_H2(
                ready,
                hap2,
                align_ont_hap2.bam
            )
        }


        // -------------------------------------------------------------------
        // Merqury
        // -------------------------------------------------------------------

        if (run_merqury) {
            merqury_hap1 = MERQURY_H1(ready, hap1, wgs_hap1_r1, wgs_hap1_r2)
            merqury_hap2 = MERQURY_H2(ready, hap2, wgs_hap2_r1, wgs_hap2_r2)
            merqury_bed_hap1 = merqury_hap1.bed
            merqury_bed_hap2 = merqury_hap2.bed
        } else {
            bed1 = params.merqury_bed_hap1 ?: params.merqury_bed
            bed2 = params.merqury_bed_hap2 ?: params.merqury_bed
            if (!bed1 || !bed2) error "ERROR: disabling Merqury requires --merqury_bed_hap1/--merqury_bed_hap2 or --merqury_bed"
            merqury_bed_hap1 = Channel.value(file(bed1, checkIfExists: true))
            merqury_bed_hap2 = Channel.value(file(bed2, checkIfExists: true))
        }

        if (!run_flagger) {
            bed1 = params.flagger_bed_hap1 ?: params.flagger_bed
            bed2 = params.flagger_bed_hap2 ?: params.flagger_bed
            if (!bed1 || !bed2) error "ERROR: disabling Flagger requires --flagger_bed_hap1/--flagger_bed_hap2 or --flagger_bed"
            flagger_hap1_bed = Channel.value(file(bed1, checkIfExists: true))
            flagger_hap2_bed = Channel.value(file(bed2, checkIfExists: true))
        } else if (run_flagger) {
            flagger_hap1_bed = flagger_hap1.bed
            flagger_hap2_bed = flagger_hap2.bed
        }


        // -------------------------------------------------------------------
        // Medaka
        //
        // Только при наличии ONT.
        // -------------------------------------------------------------------

        if (run_medaka && params.ont) {

            medaka_hap1 = MEDAKA_H1(
                ready,
                hap1,
                ont_hap1
            )

            medaka_hap2 = MEDAKA_H2(
                ready,
                hap2,
                ont_hap2
            )
        }

        ml_vcf_channels_hap1 = []
        ml_vcf_channels_hap2 = []
        if (run_deepvariant) {
            ml_vcf_channels_hap1 << deepvariant_hap1.vcf
            ml_vcf_channels_hap2 << deepvariant_hap2.vcf
        }
        if (run_pepper && params.ont) {
            ml_vcf_channels_hap1 << pepper_hap1.vcf
            ml_vcf_channels_hap2 << pepper_hap2.vcf
        }
        if (run_original_t2t) {
            ml_vcf_channels_hap1 << original_t2t_hap1.vcf
            ml_vcf_channels_hap2 << original_t2t_hap2.vcf
        }
        if (run_nextpolish2 && params.hifi) {
            ml_vcf_channels_hap1 << np2_hap1.vcf
            ml_vcf_channels_hap2 << np2_hap2.vcf
        }
        if (run_sniffles && params.hifi) {
            ml_vcf_channels_hap1 << sniffles_hifi_hap1.vcf
            ml_vcf_channels_hap2 << sniffles_hifi_hap2.vcf
        }
        if (run_sniffles && params.ont) {
            ml_vcf_channels_hap1 << sniffles_ont_hap1.vcf
            ml_vcf_channels_hap2 << sniffles_ont_hap2.vcf
        }
        if (run_cutesv && params.hifi) {
            ml_vcf_channels_hap1 << cutesv_hifi_hap1.vcf
            ml_vcf_channels_hap2 << cutesv_hifi_hap2.vcf
        }
        if (run_cutesv && params.ont) {
            ml_vcf_channels_hap1 << cutesv_ont_hap1.vcf
            ml_vcf_channels_hap2 << cutesv_ont_hap2.vcf
        }
        if (ml_vcf_channels_hap1.isEmpty()) error "ERROR: no VCF outputs available for the ML stage"

        ml_vcfs_hap1 = mix_channels(ml_vcf_channels_hap1).flatten().collect()
        ml_vcfs_hap2 = mix_channels(ml_vcf_channels_hap2).flatten().collect()
        if (workflow_mode == 'train') {
            ML_TRAIN_H1(
                ready, ml_script, ml_vcfs_hap1, truth_vcf_hap1,
                flagger_hap1_bed, merqury_bed_hap1, "hap1"
            )
            ML_TRAIN_H2(
                ready, ml_script, ml_vcfs_hap2, truth_vcf_hap2,
                flagger_hap2_bed, merqury_bed_hap2, "hap2"
            )
        } else {
            ML_PREDICT_H1(
                ready, ml_script, ml_model_hap1, ml_vcfs_hap1,
                flagger_hap1_bed, merqury_bed_hap1, "hap1"
            )
            ML_PREDICT_H2(
                ready, ml_script, ml_model_hap2, ml_vcfs_hap2,
                flagger_hap2_bed, merqury_bed_hap2, "hap2"
            )
        }


    } else {


        // ===================================================================
        // HAPLOID — оригинальная ветка (без изменений)
        // ===================================================================

        draft = Channel.value(
            file(params.draft, checkIfExists: true)
        )

        // Haploid input: bypass compression/merge for a single WGS pair;
        // merge each mate only when multiple lane files were supplied.
        wgs_lane_routes = wgs_pairs.branch {
            multiple: it[0].size() > 1
            single: it[0].size() == 1
        }
        wgs_r1_single = wgs_lane_routes.single.map { pair -> pair[0][0] }
        wgs_r2_single = wgs_lane_routes.single.map { pair -> pair[1][0] }
        wgs_r1_merged = MERGE_WGS_R1(
            wgs_lane_routes.multiple.map { pair -> tuple("wgs_R1", pair[0]) }
        )
        wgs_r2_merged = MERGE_WGS_R2(
            wgs_lane_routes.multiple.map { pair -> tuple("wgs_R2", pair[1]) }
        )
        wgs_r1 = mix_channels([wgs_r1_single, wgs_r1_merged])
        wgs_r2 = mix_channels([wgs_r2_single, wgs_r2_merged])

        if (workflow_mode == 'train') {
            truth_query = Channel.value(file(params.truth_fasta, checkIfExists: true))
            truth_vcf = TRUTH_VCF(ready, truth_script, draft, truth_query, "haploid").vcf
        }

        draft_fai = Channel.value(
            file("${params.draft}.fai", checkIfExists: true)
        )

        if (params.hifi) { hifi_reads = reads_channel(params.hifi).collect() }
        if (params.ont)  { ont_reads  = reads_channel(params.ont).collect()  }


        // -------------------------------------------------------------------
        // Alignments
        // -------------------------------------------------------------------

        if (params.hifi) {

            align_hifi = ALIGN_HIFI(
                ready,
                draft,
                hifi_reads
            )
        }

        if (params.ont) {

            align_ont = ALIGN_ONT(
                ready,
                draft,
                ont_reads
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

        if (run_deepvariant && params.hifi) {

            deepvariant = DEEPVARIANT(
                ready,
                draft,
                draft_fai,
                align_wgs.bam,
                align_wgs.bai,
                align_hifi.bam,
                true
            )

        } else if (run_deepvariant) {

            deepvariant = DEEPVARIANT(
                ready,
                draft,
                draft_fai,
                align_wgs.bam,
                align_wgs.bai,
                align_wgs.bam,
                false
            )
        }


        // -------------------------------------------------------------------
        // Pepper
        // -------------------------------------------------------------------

        if (run_pepper && params.ont) {

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

        if (run_original_t2t) {
            if (run_pepper && params.ont) {
                original_pepper = pepper.vcf
            } else {
                empty_vcf = EMPTY_VCF(deepvariant.vcf)
                original_pepper = empty_vcf.vcf
            }

            original_t2t = ORIGINAL_T2T(
                ready,
                deepvariant.vcf,
                original_pepper,
                draft,
                wgs_r1,
                wgs_r2,
                params.hifi ? hifi_reads : wgs_r1,
                params.hifi ? true : false
            )
        }


        // -------------------------------------------------------------------
        // Auto polish
        // -------------------------------------------------------------------

        if (run_auto_polish) {
            if (params.hifi) {
                auto_polish = AUTO_POLISH(
                    ready, draft, hifi_reads, wgs_r1, wgs_r2, "hifi"
                )
            } else {
                auto_polish = AUTO_POLISH(
                    ready, draft, ont_reads, wgs_r1, wgs_r2, "ont"
                )
            }
        }


        // -------------------------------------------------------------------
        // NextPolish2
        // -------------------------------------------------------------------

        if (run_nextpolish2 && params.hifi) {

            np2 = NP2(
                ready,
                draft,
                align_hifi.bam,
                wgs_r1,
                wgs_r2
            )

        }


        // -------------------------------------------------------------------
        // Sniffles
        // -------------------------------------------------------------------

        if (run_sniffles && params.hifi) {

            sniffles_hifi = SNIFFLES_HIFI(
                ready,
                draft,
                align_hifi.bam,
                align_hifi.bai,
                "hifi"
            )

        }

        if (run_sniffles && params.ont) {

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

        if (run_cutesv && params.hifi) {

            cutesv_hifi = CUTESV_HIFI(
                ready,
                draft,
                align_hifi.bam,
                align_hifi.bai,
                "hifi"
            )

        }

        if (run_cutesv && params.ont) {

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

        if (run_flagger && params.hifi) {

            flagger = FLAGGER(
                ready,
                draft,
                align_hifi.bam
            )

        } else if (run_flagger) {

            flagger = FLAGGER(
                ready,
                draft,
                align_ont.bam
            )
        }


        // -------------------------------------------------------------------
        // Merqury
        // -------------------------------------------------------------------

        if (run_merqury) {
            merqury = MERQURY(ready, draft, wgs_r1, wgs_r2)
            merqury_bed_for_ml = merqury.bed
        } else {
            if (!params.merqury_bed) error "ERROR: disabling Merqury requires --merqury_bed"
            merqury_bed_for_ml = Channel.value(file(params.merqury_bed, checkIfExists: true))
        }

        if (!run_flagger) {
            if (!params.flagger_bed) error "ERROR: disabling Flagger requires --flagger_bed"
            flagger_bed_for_ml = Channel.value(file(params.flagger_bed, checkIfExists: true))
        } else {
            flagger_bed_for_ml = flagger.bed
        }


        // -------------------------------------------------------------------
        // Medaka
        // -------------------------------------------------------------------

        if (run_medaka && params.ont) {

            medaka = MEDAKA(
                ready,
                draft,
                ont_reads
            )
        }

        ml_vcf_channels = []
        if (run_deepvariant) ml_vcf_channels << deepvariant.vcf
        if (run_pepper && params.ont) ml_vcf_channels << pepper.vcf
        if (run_original_t2t) ml_vcf_channels << original_t2t.vcf
        if (run_nextpolish2 && params.hifi) ml_vcf_channels << np2.vcf
        if (run_sniffles && params.hifi) ml_vcf_channels << sniffles_hifi.vcf
        if (run_sniffles && params.ont) ml_vcf_channels << sniffles_ont.vcf
        if (run_cutesv && params.hifi) ml_vcf_channels << cutesv_hifi.vcf
        if (run_cutesv && params.ont) ml_vcf_channels << cutesv_ont.vcf
        if (ml_vcf_channels.isEmpty()) error "ERROR: no VCF outputs available for the ML stage"

        ml_vcfs = mix_channels(ml_vcf_channels).flatten().collect()
        if (workflow_mode == 'train') {
            ML_TRAIN(
                ready, ml_script, ml_vcfs, truth_vcf,
                flagger_bed_for_ml, merqury_bed_for_ml, "haploid"
            )
        } else {
            ML_PREDICT(
                ready, ml_script, ml_model, ml_vcfs,
                flagger_bed_for_ml, merqury_bed_for_ml, "haploid"
            )
        }
    }
}
