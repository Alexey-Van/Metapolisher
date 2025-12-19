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
- **Variant calling & polishing**: Winnowmap2, bwa, Parliament2, Sniffles, TandemMapper, DeepVariant, Jasmine, Iris  
- **Annotation**: Merqury, Liftoff, RepeatMasker, Flagger  
- **Machine Learning**: CatBoost  
- **Workflow orchestration**: Nextflow, Docker/Conda  
- **QC & benchmarking**: QUAST, Merqury, hap.py  

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


