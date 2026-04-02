# Summary of copy number analysis
The data was applied on the DBS technology analysed with BLR, thus many of the parameters were adjusted accordingly <br>
If you want to use the codes for paired tumor/normal samples from the same technology, feel free to adjust the parameters.

## Steps

### Call CNVs (Titan)


Titan package was forked from [TitanCNA_10X_snakemake](https://github.com/GavinHaLab/TitanCNA_10X_snakemake) and changes on the code are summarized in [link](https://github.com/HSiga/TitanCNA_10X_snakemake/tree/working) #cbc80cc9  
Config.yaml file main modifications:
* Optional files directories were commented out.
* Gender changed male -> female
* bx_mapQual 60 -> 20
* het_minVCFQuality 100 -> 30 
* number of cores and memory for the run are adjusted
* added path for the reference genome for a faster CRAM process 

Additonal R code were cloned from [link](https://github.com/gavinha/TitanCNA/tree/master/R) and one parameter is adjusted in `haplotype.R`
* minNormQual is adjusted in file haplotype.R 100 -> 30 (in the run `haplotype2.R` was used with the modification)

The code was run on an HPC (Uppmax/Bianca) with loading pre-installed packages as in the example [srun.sh](1.TITAN/srun.sh)

Main run:
```
time snakemake -s TitanCNA.snakefile --cores n --printshellcmds >> stdout.txt 2>> stderr.txt
```
<br>

### Call CNVs (ASCAT)

The run was done using the script [cna.R](2.ASCAT/cna.R) with loading pre-installed packages as in the example [srun.sh](2.ASCAT/srun.sh).
The script was adjusted to call cram files, and to load reference SNP file from locally downloaded file.

Affymetrix SNP 6.0 (SNP6) file `GRCh38_SNP6.tsv.gz` was used as a refetence in the code for the SNPs. <br>

Versions of different packages used in the run are shown in [session.info](2.ASCAT/session_info.txt) 

Main run:
```
Rscript cna.R 1> log.txt 2> stderr.txt
```

<br>

### Call RNA-seq CNVs and compare to Titan ASCAT 

Setup conda environment specified in [cnvkit.yaml](3.cnvkit/cnvkit.yaml). 

```
micromamba env create -n cnvkit -f cnvkit.yaml
```

Analysis run using Snakemake with file [Snakefile](3.cnvkit/Snakefile). First, the file needs to be edited, providing paths to Titan, ASCAT and RNA-seq counts. 

```
micromamba activate cnvkit
snakemake --snakefile Snakefile -c 4
```

The files `out/heatmap.P19.pdf` and `out/heatmap.P18.pdf` are used in figure 2 in the article.

## Versions

Tool | Version
--- | ---
`CNVkit` | 0.9.10
