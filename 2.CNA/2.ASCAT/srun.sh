#!/bin/bash -l

#SBATCH -A sens2020007
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J CNA 

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8


source /Path/to-conda-env/bin/activate
conda activate /Path/to-conda-env

Rscript cna.R 1> log.txt 2> stderr.txt


