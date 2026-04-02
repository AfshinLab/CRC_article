#!/bin/bash -l

#SBATCH -A sens2020007
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH -J Titan

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load bioinfo-tools
module load R_packages/3.6.1

source /PATH/TO/CONDA/bin/activate 


echo -e "\n ==================================== \\n " >> log.txt
echo -e "\n ==================================== \\n " >> stderr.txt

time snakemake -s TitanCNA.snakefile --cores 16 --printshellcmds >> log.txt 2>> stderr.txt
