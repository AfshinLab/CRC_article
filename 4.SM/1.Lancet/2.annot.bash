#!/bin/bash

module load bioinfo-tools
module load samtools 
ml BioPerl/1.7.2_Perl5.26.2

# note, it worked also with the default perl in bianca
perl ~/Desktop/humam-sens2020007/mutscape_data/vvcf2maf-1.6.20/vcf2maf.pl --input-vcf P19.onlyPASS_header_fixed.vcf --output-maf P19.onlyPASS.maf --tumor-id tumor --normal-id normal --ref-fasta /castor/project/proj_nobackup/references/bwa/genome.fa --vep-path /home/humam/.conda/envs/MutScape/bin/ --ncbi-build GRCh38 --retain-fmt HPR,HPA --retain-info HPS,HPSN,HPST 

