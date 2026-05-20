# Summary of structural variants analysis

SV analysis was performed on BLR phased alignments. 

## Steps

### SV calling (LinkedSV)
Symlinked indexed phased CRAM (`*.cram` and `*.cram.crai` file) for each library into a new directory for each sample. Run LinkedSV

**N18**
```{bash}
time python /proj/sens2020007/nobackup/tools/LinkedSV/linkedsv.py --germline_mode --bam final.phased.cram --out_dir LinkedSV_N18 --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --ref_version hg38 --n_thread 16 --n_memory 5 --wgs --save_temp_files 1> stdout.log 2> stderr.log
```

**N19**
```{bash}
time python /proj/sens2020007/nobackup/tools/LinkedSV/linkedsv.py --germline_mode --bam final.phased.cram --out_dir LinkedSV_N19 --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --ref_version hg38 --n_thread 16 --n_memory 5 --wgs --save_temp_files 1> stdout.log 2> stderr.log
```

**T18**
```{bash}
time python /proj/sens2020007/nobackup/tools/LinkedSV/linkedsv.py --somatic_mode --bam final.phased.cram --out_dir LinkedSV_T18 --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --ref_version hg38 --n_thread 16 --n_memory 5 --wgs --save_temp_files 1> stdout.log 2> stderr.log 
```

**T19**
```{bash}
time python /proj/sens2020007/nobackup/tools/LinkedSV/linkedsv.py --somatic_mode --bam final.phased.cram --out_dir LinkedSV_T19 --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --ref_version hg38 --n_thread 16 --n_memory 5 --wgs --save_temp_files 1> stdout.log 2> stderr.log 
```

LinkedSV code modified to accept custom parameter "--n_memory"  that allow for higher memory usage when sorting per barcode.   

### SV calling (NAIBR)

SV calling with NAIBR was performed as part of BLR pipeline.

### Combine LinkedSV calls per patient

Symlink small deletions (file *.small_deletions.bedpe) and filtered large SVs (file *.filtered_large_svcalls.bedpe) into new directory

```{bash}
base="" # Path to LinkedSV output
for i in N18 T18 N19 T19
do
	echo $i
	ln -s "${base}/$i/LinkedSV_${i}/final.phased.cram.filtered_large_svcalls.bedpe" "${i}.filtered_large_svcalls.b
	ln -s "${base}/$i/LinkedSV_${i}/final.phased.cram.small_deletions.bedpe" "${i}.small_deletions.bedpe"
done
```

Concatenate calls per patient

```{bash}
for i in 18 19
do
    echo $i
    cat N$i*.bedpe T$i*.bedpe | grep -v "^#" | grep "^chr" | cut -f -7 | sort -k 1,1 -k2,2n -k4,4 -k5,5n | uniq > P$i.linkedsv.candidates.no_type.bedpe
done
```

Use Python script [add_type.py](add_type.py) to add breakpoint direction information for NAIBR

```{bash}
cat P19.linkedsv.candidates.no_type.bedpe | python add_type.py > P19.linkedsv.candidates.bedpe
cat P18.linkedsv.candidates.no_type.bedpe | python add_type.py > P18.linkedsv.candidates.bedpe
```

### Combine NAIBR calls per patient

Symlink NAIBR BEDPEs to the current directory and add a sample prefix (e.g. N18, N19, T18, T19). Concatenate calls per patient

```{bash}
cat ?18.naibr_sv_calls.bedpe | grep -v "^#" | grep PASS | cut -f 1-6,12 | tr ";" "\t" | tr "=" "\t" | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$18,$14}' > P18.naibr.candidates.bedpe 
cat ?19.naibr_sv_calls.bedpe | grep -v "^#" | grep PASS | cut -f 1-6,12 | tr ";" "\t" | tr "=" "\t" | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$18,$14}' > P19.naibr.candidates.bedpe
```

### Combine all candidate calls per patient
Concatenate all calls from LinkedSV and NAIBR per patient, sort and filter out duplicates

```{bash}
cat P18.*.candidates.bedpe | sort -k 1,1 -k2,2n -k4,4 -k5,5n | uniq > P18.candidates.bedpe
cat P19.*.candidates.bedpe | sort -k 1,1 -k2,2n -k4,4 -k5,5n | uniq > P19.candidates.bedpe
```

### Re-score calls with NAIBR
For each sample, the patient candidate calls were re-scored using NAIBR

**Configurations**
Example for N18.

Update arguments for other samples:
- `bam_file`: path to phased CRAM for sample
- `candidates`: path to candidate calls for patient
- `prefix`: output prefix for sample
- `outdir`: output directory, e.g. sample name

```
# minimum mapping quality (default=40)
min_mapq=40

# input bam file
bam_file=/proj/nobackup/sens2020007/analysis_nobackup/RERUNS_WITH_CHRX/N18/final.phased.cram

# prefix
prefix=N18.NAIBR_SVs

# output directory (default=.)
outdir=N18

# list of intervals not to be intcluded in analysis (default=None)
#blacklist=

# list in BEDPE format of novel adjacencies to be scored by NAIBR (default=None)
candidates=create_candidates/P18.candidates.bedpe

# maximum distance between read-pairs in a linked-read (default=10000)
d=10000

# minimum size of structural variant (default=2*lmax)
min_sv=2000

# number of threads (default=1)
threads=16

# minimum number of barcode overlaps supporting a candidate NA (default = 3)
k=3

# minimum length of linked-read fragment to consider (default=2*lmax)
#min_len=

# minimum nr of reads in linked-read fragment for it to be considered (default=2)
min_reads=2

# minimum number of discordant reads required (default = 2)
min_discs=2
```

**Run NAIBR**
Example for N18

```{bash}
naibr N18.naibr.config
```

### Prepare SV blacklists

Downloaded 10x blacklists from: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/references

```{bash}
mkdir sv_blacklists
cd sv_blacklists
wget -O GRCh38_10x_sv_blacklist.bed https://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed
wget -O GRCh38_10x_segdups.bedpe https://cf.10xgenomics.com/supp/genome/GRCh38/segdups.bedpe
```

Get chromosome lengths. The blacklist is from 10X Genomics LongRanger pipeline and contains `alt` contigs that we need to remove.

```{bash}
wget -O GRCh38.genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
```

Get chromosome names

```{bash}
cut -f 1 GRCh38.genome > GRCh38.chroms.list
```

Select only contigs from GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set

```{bash}
grep -w -f GRCh38.chroms.list GRCh38_10x_sv_blacklist.bed > GRCh38_10x_sv_blacklist.select.bed 
```

Expand BED 10kb using BEDtools

```{bash}
bedtools slop -i GRCh38_10x_sv_blacklist.select.bed -g GRCh38.genome -b 10000 > GRCh38_10x_sv_blacklist.slop_10kbp.bed
```


### Filter re-scored SVs
Symlink NAIBR SVs

```{bash}
base="" #Path to dir with NAIBR re-scored SVs
for i in N18 T18 N19 T19
do
    echo $i
    ln -s "${base}/${i}/${i}.NAIBR_SVs.vcf" "${i}.vcf"
done
```

Merge VCF per patient

```{bash}
ref="/proj/nobackup/sens2020007/references/bwa/genome.fa" # Path to genome reference

for p in 18 19
do
	
	bcftools reheader -s <(echo "T${p}") "T${p}.vcf" | bcftools sort -o "T${p}.vcf.gz"
	tabix "T${p}.vcf.gz"
        
	bcftools reheader -s <(echo "N${p}") "N${p}.vcf" | bcftools sort -o "N${p}.vcf.gz"
        tabix "N${p}.vcf.gz"
	
	bcftools merge -m none -f PASS "T${p}.vcf.gz" "N${p}.vcf.gz" | bgzip -c > "P${p}.merge.vcf.gz"

	tabix "P${p}.merge.vcf.gz"

	truvari collapse -i "P${p}.merge.vcf.gz" -o "P${p}.truvari_merge.vcf" -c "P${p}.truvari_collapsed.vcf" -f $ref \
	 --chain -k maxqual --pctsim 0 --pctovl 0.8 --pctsize 0.8 --refdist 1000 -S 100_000 --passonly

done
```

Filter calls for size and against blacklists

```{bash}
blacklist_bed="sv_blacklists/GRCh38_10x_sv_blacklist.slop_10kbp.bed"
segdups="sv_blacklists/GRCh38_10x_segdups.bedpe"

for p in 18 19
do
	# Filter blacklist and size 2-100kb
	SURVIVOR filter "P${p}.truvari_merge.vcf" $blacklist_bed 2000 100000 -1 -1 "P${p}.truvari_merge.filt.tmp.vcf" > "P${p}.truvari_merge.filt.vcf.log"
	
	# Filter segdups
	lsvtool intersect_bedpe "P${p}.truvari_merge.filt.tmp.vcf" $segdups -o "P${p}.truvari_merge.filt.vcf" -d 20_000 2>> "P${p}.truvari_merge.filt.vcf.log"

done
```

Get somatic SV calls

```{bash}
for p in 18 19
do
	bcftools sort "P${p}.truvari_merge.filt.vcf" -o "P${p}.truvari_merge.filt.vcf.gz"

	tabix "P${p}.truvari_merge.filt.vcf.gz"

	bcftools +setGT "P${p}.truvari_merge.filt.vcf.gz" -- -t . -n 0 | bcftools view -i 'GT[1]="ref"' -o "P${p}.truvari_merge.filt.somatic.vcf"

done
```

<br>

## Haploplot visualization
For visulaising the selected SVs [haploplot](https://github.com/HSiga/haploplot) was used, with the following parameters:

Example of the running code:

```{bash}
python haploplot.py -p chr5:19320000-19480000  -s1 ../N18.cram -l1 N18 -s2 ../T18.cram -l2 T18 -c1 FFA382 -w 5000 -r chr5:19384245-19419473 -f 30000 -o P18_del_chr5_ann.png -t del -gff3 sorted.gff3.gz -ref genome.fa 1> sv_del_chr5.log 2>&1 &

python haploplot.py -p chr4:95183000-95352000 -s1 ../N19.cram -l1 N19 -s2 ../T19.cram -l2 T19 -c1 FFA382 -w 5000 -r chr4:95225267-95310030 -f 30000 -o P19_del_chr4_ann.png -t del  -gff3 sorted.gff3.gz -ref genome.fa 1> sv_del_chr4.log 2>&1 &

python haploplot.py -p chr20:38393000-38456000 -s1 ../N19.cram -l1 N19 -s2 ../T19.cram -l2 T19 -c1 FFA382 -w 500 -r chr20:38409153-38440198 -f 500 -o P19_inv_chr20_ann.png -t inv  -gff3 sorted.gff3.gz -ref genome.fa 1> sv_inv_chr20.log 2>&1

wait

```

<br>

## Sashimi plot
Sashimi plot used in the supplementary was generated using the code [sashimi-plot.py](sashimi-plot.py). 

Example of the running code:

```{bash}
python sashimi.py --bam1 N18.cram --bam2 T18.cram --l1 N18 --l2 T18 --region chr5:19370000-19435000 --highlight chr5:19384245-19419473,chr5:19375138-19376301 --highlight_color fff6d9 --output sash_P18_del.png --swap_hp bam1 
```
  
<br>

## Versions

Tool | Link | Version
--- | --- | ---
`NAIBR`| https://github.com/pontushojer/NAIBR | v0.5.1
`LinkedSV` | https://github.com/HSiga/LinkedSV | commit #84186a9 
`haploplot` | https://github.com/HSiga/haploplot/ | 0.1.0 
`BCFtools` | | v1.17
`BEDtools`| |  v2.30.0
`Truvari` | | v3.2.0
`SURVIVOR` | | 1.0.7
