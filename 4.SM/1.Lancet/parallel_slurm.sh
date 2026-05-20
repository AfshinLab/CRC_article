#!/bin/bash -l

#SBATCH -A sens2020007
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J lancet_chrs_P18
#SBATCH --array=1-22

echo "started lancet ${SLURM_ARRAY_TASK_ID}"

time lancet --linked-reads --primary-alignment-only --tumor T18_inputs/chr${SLURM_ARRAY_TASK_ID}.calling.phased.bam --normal N18_inputs/chr${SLURM_ARRAY_TASK_ID}.calling.phased.bam --ref /proj/sens2020007/nobackup/references/bwa/genome.fa --reg "chr${SLURM_ARRAY_TASK_ID}" --num-threads 16 > runs_vcf/chr${SLURM_ARRAY_TASK_ID}_out.vcf 2> runs_txt/log_chr${SLURM_ARRAY_TASK_ID}_lancet.txt

echo "lancet is done! ${SLURM_ARRAY_TASK_ID}"
