# CRC_article
Collection of scripts for analyzing DBS data from colorectal cancer samples

### Introduction
Paired tumor/normal whole-genome samples were sequenced with linked-read sequencing technology (DBS technology). 

In this repo we will summarize the analysis steps with the main scripts, versions and parameters applied in every section.

 
### Main steps of the analysis

Reference genome (GRCh38): https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

1. [BLR analysis](1.BLR)
2. [Copy number analysis](2.CNA)
    - TitanCNA (linked-reads workflow)
    - ASCAT (maftools workflow)
    - Cnvkit
3. [Structural variant analysis](3.SV)
    - NAIBR
    - LinkedSV
4. [Small (somatic) mutation analysis](4.SM)
    - Lancet
    - Maftools
5. [RNAseq analysis](5.RNA)
    - R scripts


