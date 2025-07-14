# Single-cell immune cell atlas with multi-layer omics data in Japanese

## Overview
We constructed multi-omics immune cell atlas, <ins>**O**</ins>saka <ins>**A**</ins>tla<ins>**s**</ins> of <ins>**I**</ins>mmune Cell<ins>**s**</ins> (OASIS), from 235 Japanese including COVID-19 patients and healthy subjects. The 5’ single-cell transcriptomics data profiled over 1,500,000 peripheral blood mononuclear cells. OASIS links these single-cell transcriptomics data with host genetics, plasma proteomics, and metagenomics data.\
\
This repository provides code used for single-cell eQTL mapping, colocalization analysis, PRS analysis, megatenomics analysis, and single-cell deconvolution of various somatic mutations.\
The code used for the HLA and genome-wide association analysis of immune repertoires is available at [https://github.com/tatsuhikonaito/OASIS_HLATCR](https://github.com/tatsuhikonaito/OASIS_HLATCR).

![OverView](https://github.com/user-attachments/assets/92d898a9-1bdd-4643-935f-dda6fd3b4a72)

Edahiro R, Sato G, Naito T, et al, "Deciphering state-dependent immune features from multi-layer omics data at single-cell resolution", Nature Genetics in press. 

## Data availability
Raw sequencing data of scRNA-seq and the protein expression matrix are available at the Japanese Genotype-phenotype Archive (JGA) with accession codes [JGAS000783](https://ddbj.nig.ac.jp/search/entry/jga-study/JGAS000783)/[JGAD000925](https://ddbj.nig.ac.jp/search/entry/jga-dataset/JGAD000925). A part of the raw scRNA-seq data (nCOVID-19 = 73, nControl = 75) has already been deposited and is available under controlled access at JGA with accession codes JGAS000593/JGAS000543/JGAD000662/JGAD000722. All the raw sequencing data of scRNA-seq can also be accessed through application at the NBDC with the accession code [hum0197](https://humandbs.biosciencedbc.jp/en/hum0197-latest). Genotype data of the subjects are available at European Genome-Phenome Archive (EGA) with the accession code [EGAS00001008016](https://ega-archive.org/studies/EGAS00001008016). The sc-eQTL summary statistics are available at NBDC accession ID [hum0197](https://humandbs.biosciencedbc.jp/en/hum0197-latest), and also available at an interactive browser: [JOB (<ins>**J**</ins>apan <ins>**O**</ins>mics <ins>**B**</ins>rowser)](https://japan-omics.jp/).

## Contact
Ryuya Edahiro: r.edahiro_at_imed3.med.osaka-u.ac.jp, redahiro_at_sg.med.osaka-u.ac.jp 
