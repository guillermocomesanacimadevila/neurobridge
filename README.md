<p align="center">
  <img src="https://github.com/user-attachments/assets/7ff124df-c97d-44e2-b51a-6617a80ad61d"
       alt="nf-core-neurobridge_logo_light"
       width="1600">
</p>

---

## Authors

<p align="center">

<strong>Guillermo Comesaña Cimadevila<sup>1,2</sup></strong><br>
Dervis Salih<sup>1,3</sup> · Nicholas Bray<sup>2</sup> · Emily Simmonds<sup>1</sup> · Valentina Escott-Price<sup>1,2</sup>

<br><br>

<sup>1</sup> UK Dementia Research Institute, UK 
<sup>2</sup> MRC Centre for Neuropsychiatric Genetics & Genomics, Cardiff University, Cardiff, UK  
<sup>3</sup> Department of Neuroscience, Physiology and Pharmacology, University College London (UCL), London, UK  

</p>

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

### Build Docker image

```bash
cd neurobridge/
docker build -t neurobridge:1 -f env/Dockerfile env/
```
