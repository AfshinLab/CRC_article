# Summary of copy number analysis

## Steps

### Call CNVs (Titan)

TODO

### Call CNVs (ASCAT)

TODO

### Call RNA-seq CNVs and compare to Titan ASCAT 

Setup conda environment specified in `cnvkit.yaml`. 

```
conda env create -n cnvkit -f cnvkit.yaml
```

Analysis run using Snakemake with file `Snakefile`. First, the file needs to be edited, providing paths to Titan, ASCAT and RNA-seq counts. 

```
conda activate cnvkit
snakemake
```

The files `out/heatmap.P19.pdf` and `out/heatmap.P18.pdf` are used in figure 2 in the article.

## Versions

Tool | Link | Version
--- | --- | ---
`CNVkit` | | 0.9.10
