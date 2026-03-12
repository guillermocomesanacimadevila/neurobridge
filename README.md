<p align="center">
  <img src="https://github.com/user-attachments/assets/7ff124df-c97d-44e2-b51a-6617a80ad61d"
       alt="nf-core-neurobridge_logo_light"
       width="1600">
</p>

---

## Introduction

**Neurobridge** is a two-stage bioinformatics pipeline which takes two genome-wide association study (GWAS) summary statistics corresponding to two distinct traits (typically trait pairs where a pleiotropic relationship may be hypothesised). The pipeline spans analyses from estimating global heritability to functional analyses using bulk and single-cell xQTL datasets. In **Stage I**, neurobridge performs GWAS QC and computes SNP heritability, polygenicity and discoverability, global genetic correlation (across multiple methods), and local genetic correlation. It then identifies loci showing evidence of pleitropic association and performs Bayesian colocalisation, and fine-mapping to prioritise candidate causal variants. Following this stage, loci can be annotated using **FUMA** (performed manually) to map variants to genes and obtain basic functional annotation. In **Stage II**, neurobridge focuses on functional interpretation of the prioritised loci by integrating bulk and single-cell xQTL datasets and performing gene prioritisation (via SMR). In essence, Stage II adds multi-layer biological evidence to the genetic findings from Stage I.

---

## Authors

**Guillermo Comesaña Cimadevila¹ ²** · Dervis Salih¹ ³ · Nicholas Bray² · Emily Simmonds¹ · Valentina Escott-Price¹ ²

¹ UK Dementia Research Institute, UK  
² MRC Centre for Neuropsychiatric Genetics & Genomics, Cardiff University, Cardiff, UK  
³ Department of Neuroscience, Physiology and Pharmacology, University College London, London, UK  

---

[![Python](https://img.shields.io/badge/Python-3.12%2B-3776AB)]() 
[![R](https://img.shields.io/badge/R-4.4%2B-276DC3)]() 
[![Run with Conda](https://img.shields.io/badge/Run%20with-Conda-44A833)]() 
[![Nextflow DSL2](https://img.shields.io/badge/Nextflow-DSL2-23aa62)](https://www.nextflow.io/) 
[![Docker](https://img.shields.io/badge/Container-Docker-2496ED)](https://www.docker.com/) 
[![Reference](https://img.shields.io/badge/Docs-auto--generated-green)](zenodo)

---

## Pipeline

<p align="center">
  <img src="https://github.com/user-attachments/assets/a15c8413-a93c-43b7-b647-b0fa5a38def2"
       alt="pipeline"
       width="1200">
</p>

---

## Steps

### Download reference data from Zenodo

> Download **[here](https://doi.org/10.5281/zenodo.18986935)** 

---

### Build Docker image

```bash
cd neurobridge/
```

```bash
docker build -t neurobridge:1 -f env/Dockerfile env/
```

---

### Run neurobridge! (with Docker)

```bash
nextflow run main.nf \
  -profile docker \
  -c conf/local/nextflow.config \
  --input assets/gwas.tsv \
  --pairs assets/ldsc_pairs.tsv \
  --outdir results
```

# Packages Used

## Heritability & Genetic Correlation
- **[LDSC](https://github.com/bulik/ldsc)** – Linkage Disequilibrium Score Regression  
- **[SumHer](https://github.com/zietzm/sumher_rs)** – Summary-based Heritability Estimation  
- **[HDL / HDL-L](https://github.com/zhenin/HDL)** – High-Definition Likelihood (Global and Local)  
- **[MiXeR](https://github.com/precimed/mixer)** – Mixture of Regressions  
- **[LAVA](https://github.com/josefin-werme/LAVA)** – Local Analysis of [co]Variant Association  

## Cross-trait Enrichment
- **[condFDR / conjFDR](https://github.com/alexploner/cfdr.pleio)** – Conditional / Conjunctional False Discovery Rate  

## GWAS Processing & Gene-Level Analysis
- **[PLINK 2.0](https://www.cog-genomics.org/plink/2.0/)** – Whole-genome Association Analysis Toolset  
- **[MAGMA](https://ctg.cncr.nl/software/magma)** – Multi-marker Analysis of GenoMic Annotation  

## Fine-mapping & Colocalisation
- **[COLOC](https://cran.r-project.org/web/packages/coloc/)** – Bayesian Colocalisation Analysis  
- **[SuSiE](https://github.com/stephenslab/susieR)** – Sum of Single Effects

## Functional Annotation & QTL Integration
- **[FUMA](https://fuma.ctglab.nl/)** – Functional Mapping and Annotation  
- **[SMR](https://yanglab.westlake.edu.cn/software/smr/)** – Summary-data-based Mendelian Randomization (GWAS–QTL integration framework)  

---

## Get Started!

1. Install [`Nextflow`](https://www.nextflow.io/) `(>=25.10.0)`

2. Install [`Docker`](https://www.docker.com/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/)

3. Clone repo and run!
   
```bash
git clone https://github.com/guillermocomesanacimadevila/neurobridge.git
```

```bash
cd neurobridge/
```

4. Run end-to-end!

```bash
nextflow run main.nf \
  -profile docker \
  -c conf/local/nextflow.config \
  --input assets/gwas.tsv \
  --pairs assets/ldsc_pairs.tsv \
  --outdir results
```

---

** Note ** 

> You can also run each method within [`workflows/neurobridge`], individually like so =>

```bash
nextflow run workflows/neurobridge/main_<input_method>.nf \
  -profile docker \
  -c conf/local/nextflow.config \
  --input assets/gwas.tsv \
  --pairs assets/ldsc_pairs.tsv \
  --outdir results -resume
```
