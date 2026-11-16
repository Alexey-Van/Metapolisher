# MetaPolisher

## Goal
To create and experimentally test a system for automated polishing of genome assemblies **MetaPolisher** based on CatBoostClassifier.

## Relevance
Existing automatic pipelines (Racon, Medaka, DeepPolisher, etc.) do not always correctly detect and fix structural variants (SV/SNV). This leads to false positive results and reduces the reliability of the final assembly.

**MetaPolisher** combines modern approaches to polishing, variant calling, genome annotation, and machine learning to obtain a genome with quality comparable to CHM13 while reducing human labor hours.

The development of such an approach will allow:
- To increase the accessibility and convenience of polishing.
- To accelerate the analysis process.
- To reduce the proportion of false positive changes in genome assemblies.  
This is especially relevant for **T2T** and **Pangenome** projects.

---

## Project Tasks
- Perform classification and annotation of the genome.  
- Collect and integrate a set of tools for variant detection and assembly error evaluation.  
- Combine datasets into a single table and perform classification using **CatBoost** (gradient boosting on decision trees).  
- Create a prototype of an automated pipeline in a reproducible environment (**Nextflow**).  
- Evaluate the effectiveness of the approach on benchmark data (**HG002, CHM13**) using metrics:
  - QV (Quality Value)  
  - k‑mer completeness  
  - false correction rate  
- Compare the results of the system with existing polishing methods.  

---

## Technologies and Tools
- **Variant calling & polishing**: Winnowmap, bwa-mem3, Parliament2, Sniffles, TandemMapper, DeepVariant, Jasmine, Iris  
- **Annotation**: Merqury, Liftoff, RepeatMasker, Flagger  
- **Machine Learning**: CatBoost  
- **Workflow orchestration**: Nextflow, Docker/Conda  
- **QC & benchmarking**: QUAST, Merqury, hap.py  

---

## Nextflow workflow

This section covers setup, input preparation, prediction, training, configuration, and output files. Run commands from the repository root so the relative paths below resolve as written.

### 1. Requirements and container setup

Install Nextflow with DSL2 support and Docker, and make sure the Docker daemon is running. The local profile enables Docker. Build the workflow images once, and rebuild them after changing a Dockerfile:

```bash
bash scripts/build_container.sh
```

The build script also pulls the external DeepVariant, Pepper, and Flagger images, so it needs network access. The workflow checks for the images required by the selected tools before starting.

The local profile is configured for a large machine (up to 32 CPUs and 300 GB RAM per task). Adjust the `local` profile resource defaults in `scripts/nextflow.config` for a smaller host before running. The `slurm` profile selects the Slurm executor and tool-specific resources. It requires Docker and the built images to be accessible on compute nodes; clusters that require Apptainer/Singularity need site-specific container configuration.

### 2. Prepare input files

Every run requires paired WGS reads and at least one long-read input. Supply HiFi reads, ONT reads, or both:

| Inputs | Workflow routes |
| --- | --- |
| WGS + ONT | WGS alignment/calling plus ONT alignment and ONT-compatible callers/polishers |
| WGS + HiFi | WGS + HiFi alignment, hybrid DeepVariant, and HiFi-compatible tools |
| WGS + HiFi + ONT | Both long-read routes; callers use alignments for their matching platform |

`wgs`, `hifi`, and `ont` accept a FASTQ file, a glob, or a directory. WGS must be paired-end. File names must end in `_R1`/`_R2` or `_1`/`_2` before the FASTQ suffix; lane suffixes such as `_001` are supported. R1 and R2 must have matching sample and lane names. Multiple WGS lanes are merged before alignment; long-read files are passed together to their aligner. Provide data for one biological sample per run.

Examples of supported names:

```text
sample_R1_001.fastq.gz    sample_R2_001.fastq.gz
sample_R1_002.fastq.gz    sample_R2_002.fastq.gz
sample_1.fq.gz            sample_2.fq.gz
```

Create a FASTA index for the draft and, in diploid mode, both haplotype assemblies. The `.fai` must be beside the FASTA and have the same path plus `.fai`:

```bash
samtools faidx /data/draft.fa
samtools faidx /data/hap1.fa
samtools faidx /data/hap2.fa
```

Training also needs a truth assembly. Use the draft as REF and the known correct assembly as QRY. SyRI works best when corresponding chromosomes have matching contig names, so make the reference and truth FASTA headers consistent. The truth FASTA does not need a `.fai`; the truth stage creates a task-local index for its staged REF.

### 3. Configure the run

Create a JSON params file. Use absolute paths for input files and models; quote paths and use JSON booleans (`true`/`false`). `wgs` is preferred over the legacy aliases `illumina` and `mgi`.

Common parameters:

| Parameter | Purpose |
| --- | --- |
| `workflow_mode` | `predict` or `train`; defaults to `predict` |
| `ploidy` | `haploid` (default) or `diploid` |
| `draft` | Haploid draft assembly; used as the REF for calling and training |
| `hap1`, `hap2` | Diploid draft haplotype assemblies; use these instead of `draft` |
| `wgs` | Paired WGS FASTQ file, glob, or directory; required in every mode |
| `hifi`, `ont` | HiFi and/or ONT FASTQ file, glob, or directory; at least one is required |
| `outdir` | Results directory; defaults to `results` |
| `ml_threshold` | Prediction probability threshold from 0 to 1; defaults to `0.5` |

Prediction model parameters:

| Ploidy | Parameters |
| --- | --- |
| Haploid | `ml_model`: path to the required CatBoost `.cbm` file |
| Diploid | Use one shared `ml_model`, or provide both `ml_model_hap1` and `ml_model_hap2` |

Training truth parameters:

| Ploidy | Parameters |
| --- | --- |
| Haploid | `truth_fasta`: known correct assembly corresponding to `draft` |
| Diploid | `truth_fasta_hap1` and `truth_fasta_hap2`, corresponding to `hap1` and `hap2` |

### 4. Run prediction

Example haploid prediction with WGS, HiFi, and ONT:

```json
{
  "workflow_mode": "predict",
  "ploidy": "haploid",
  "draft": "/data/draft.fa",
  "wgs": "/data/wgs/*.fastq.gz",
  "hifi": "/data/hifi/*.fastq.gz",
  "ont": "/data/ont/*.fastq.gz",
  "ml_model": "/models/metapolisher.cbm",
  "ml_threshold": 0.5,
  "outdir": "results"
}
```

For HiFi-only or ONT-only runs, omit the other long-read parameter. To run the workflow:

```bash
nextflow run scripts/MetaPolisher.nf \
  -c scripts/nextflow.config \
  -profile local \
  -params-file params.json
```

Use `-resume` after a failed or interrupted run to reuse completed Nextflow tasks:

```bash
nextflow run scripts/MetaPolisher.nf \
  -c scripts/nextflow.config \
  -profile local \
  -params-file params.json \
  -resume
```

For a Slurm environment, replace `-profile local` with `-profile slurm` and add `-with-docker` if Docker is the runtime available on your compute nodes. The current container check uses the Docker CLI and expects the configured images to be available to compute nodes. Clusters that prohibit Docker need a site-specific Apptainer/Singularity setup and a matching container check.

Example diploid prediction with different models for each haplotype:

```json
{
  "workflow_mode": "predict",
  "ploidy": "diploid",
  "hap1": "/data/hap1.fa",
  "hap2": "/data/hap2.fa",
  "wgs": "/data/wgs/*.fastq.gz",
  "hifi": "/data/hifi/*.fastq.gz",
  "ml_model_hap1": "/models/hap1.cbm",
  "ml_model_hap2": "/models/hap2.cbm",
  "outdir": "results"
}
```

In diploid mode, the workflow assigns WGS and long reads to haplotypes before alignment and runs the final ML stage separately on hap1 and hap2. HiFi is preferred by tools that accept either long-read type; NextPolish2 runs only when HiFi is supplied, while Pepper and Medaka run only when ONT is supplied.

### 5. Train a model

Training uses the same read and draft inputs as prediction, plus the known correct truth assembly. The workflow runs `my_scripts/truth_paftools_syri.sh` to align the truth FASTA to the draft REF with minimap2, SyRI, and paftools, then trains CatBoost from the enabled VCF callers and Flagger/Merqury annotations. Do not set `ml_model` for training.

Example haploid training params:

```json
{
  "workflow_mode": "train",
  "ploidy": "haploid",
  "draft": "/data/draft.fa",
  "truth_fasta": "/data/known_correct.fa",
  "wgs": "/data/wgs/*.fastq.gz",
  "ont": "/data/ont/*.fastq.gz",
  "outdir": "training_results"
}
```

Run it with the same Nextflow command, using the training params file. For diploid training, provide separate truth FASTAs for each haplotype:

```json
{
  "workflow_mode": "train",
  "ploidy": "diploid",
  "hap1": "/data/hap1.fa",
  "hap2": "/data/hap2.fa",
  "truth_fasta_hap1": "/data/truth_hap1.fa",
  "truth_fasta_hap2": "/data/truth_hap2.fa",
  "wgs": "/data/wgs/*.fastq.gz",
  "hifi": "/data/hifi/*.fastq.gz",
  "outdir": "training_results"
}
```

Diploid training writes one model for each haplotype. Use those files later with `ml_model_hap1` and `ml_model_hap2` in diploid prediction. The workflow needs at least one enabled VCF-producing caller. For training, each class must have at least three candidate variants so the script can create stratified train, validation, and test splits.

### 6. Select tools and supply precomputed annotations

All tool flags default to `true`. Disable tools independently in the params JSON when they are not needed:

| Parameter | Tool/stage | Conditions or notes |
| --- | --- | --- |
| `run_deepvariant` | DeepVariant VCF calling | Uses WGS; with HiFi it uses hybrid calling |
| `run_pepper` | PEPPER VCF calling | Runs only with ONT |
| `run_original_t2t` | Original T2T VCF pipeline | Requires DeepVariant |
| `run_auto_polish` | Automated polishing | Uses HiFi when available, otherwise ONT |
| `run_nextpolish2` | NextPolish2 | Runs only with HiFi |
| `run_sniffles` | Sniffles SV calling | Runs separately for each supplied long-read platform |
| `run_cutesv` | cuteSV calling | Runs separately for each supplied long-read platform |
| `run_flagger` | Flagger BED annotation | If disabled, provide a precomputed BED |
| `run_merqury` | Merqury BED annotation | If disabled, provide a precomputed BED |
| `run_medaka` | Medaka polishing | Runs only with ONT |

Example: disable PEPPER, cuteSV, and Medaka:

```json
{
  "run_pepper": false,
  "run_cutesv": false,
  "run_medaka": false
}
```

The same booleans can be passed on the Nextflow command line, for example `--run_pepper false --run_cutesv false --run_medaka false`.

At least one applicable VCF-producing stage must remain enabled: DeepVariant, PEPPER with ONT, Original T2T, NextPolish2 with HiFi, Sniffles, or cuteSV. `run_auto_polish` and `run_medaka` produce polishing outputs but do not add VCFs to the final ML feature set. If Flagger or Merqury is disabled, provide its BED because both tracks are used by the ML script:

```json
{
  "run_flagger": false,
  "flagger_bed": "/data/annotations/flagger.bed",
  "run_merqury": false,
  "merqury_bed": "/data/annotations/merqury.bed"
}
```

For diploid mode, use `flagger_bed_hap1`, `flagger_bed_hap2`, `merqury_bed_hap1`, and `merqury_bed_hap2` for separate tracks. A shared `flagger_bed` or `merqury_bed` can be used for both haplotypes. `run_original_t2t` cannot be enabled when `run_deepvariant` is false.

Precomputed BED files must use contig names and coordinates for the matching draft or haplotype REF. Keep the VCF-producing tools and BED annotation sources consistent between training and prediction where possible, so the model receives comparable feature columns.

You can also override params on the command line, for example `--outdir results_run2` or `--ml_threshold 0.7`. The JSON file is more convenient for recording a complete run configuration.

### 7. Find and interpret results

The workflow writes under `outdir` (default `results/`):

| Output | Location |
| --- | --- |
| Prediction VCF and full variant feature table | `ml/<prefix>/<prefix>.ml.vcf` and `<prefix>.ml.variants.tsv` |
| Trained model, variant table, hyperparameter report | `ml_train/<prefix>/model_<prefix>.cbm`, `<prefix>.variants.tsv`, and `<prefix>_hyperparams.tsv` |
| Generated truth VCF and index | `truth_vcf/<prefix>/truth.vcf.gz` and `truth.vcf.gz.csi` |
| Intermediate caller, alignment, and annotation outputs | Tool-specific folders such as `deepvariant/`, `align/`, `flagger/`, and `merqury/` |
| Nextflow timeline, report, and trace | `timeline.html`, `report.html`, and `trace.txt` |

For example, haploid training with `outdir: "training_results"` creates `training_results/ml_train/haploid/model_haploid.cbm`. Diploid training creates `training_results/ml_train/hap1/model_hap1.cbm` and `training_results/ml_train/hap2/model_hap2.cbm`.

The prediction VCF contains probability-filtered SNVs; the variants TSV contains the full feature table used by ML. In training output, `model_<prefix>.cbm` is the CatBoost model. Keep the corresponding params file, model, and hyperparameter TSV together to make later prediction runs reproducible.

The container recipes pin the main tool versions: bwa-mem3 0.17.0, samtools/bcftools 1.25, minimap2 2.31, Winnowmap 2.03, Racon 1.5.0, Sniffles 2.8.1, cuteSV 2.1.4, Merqury 1.4.1, Merfin 1.1, NextPolish2 0.2.2, Medaka 2.2.2, SyRI 1.8.2, and k8 1.2. The external images are DeepVariant 1.10.0 and Flagger v1.2.1; Pepper uses r0.8.

---

## Quality Metrics
- **QV (Quality Value)** — assembly accuracy in the Phred scale  
- **k‑mer completeness** — completeness across the k‑mer spectrum  
- **False correction rate** — proportion of false corrections  

---

## Progress Checklist

- [x] Genome classification and annotation  
- [x] Integration of tools for variant calling and error evaluation  
- [x] Combining data into a single table  
- [x] Implementation of classification with CatBoost  
- [ ] Creation of a pipeline prototype in Nextflow  
- [ ] Testing on benchmark data (HG002, CHM13)  
- [ ] Calculation of quality metrics (QV, k‑mer completeness, false correction rate)  
- [ ] Comparison with existing polishing methods  
- [ ] Preparation of the final report and visualizations  
 

---


