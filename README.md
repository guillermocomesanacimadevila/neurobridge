<p align="center">
  <img src="https://github.com/user-attachments/assets/7ff124df-c97d-44e2-b51a-6617a80ad61d"
       alt="nf-core-neurobridge_logo_light"
       width="1600">
</p>

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
  <img src="https://github.com/user-attachments/assets/318d1e31-0309-4d91-9d9b-c4b837a6e963"
       alt="pipeline"
       width="1200">
</p>

---

## Steps

### Download reference data 

=> Zenodo link...

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
nextflow run main.nf -profile docker \
  --input assets/gwas.tsv \
  --pairs assets/ldsc_pairs.tsv \
  --outdir results \
  -resume
```

