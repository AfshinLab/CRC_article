#!/bin/bash -l

#SBATCH -A sens2020007
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J lancet_chrX

echo "started lancet X "

time lancet --linked-reads --primary-alignment-only --tumor T18_inputs/chrX.calling.phased.bam --normal N18_inputs/chrX.calling.phased.bam --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --reg "chrX" --num-threads 16 > runs_vcf/chrX_out.vcf 2> runs_txt/log_chrX_lancet.txt

echo "lancet is done! X"
