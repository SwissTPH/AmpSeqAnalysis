#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/scicore/home/pothin/golmon00/JOB_OUT/qc_%A_%a.log
#SBATCH --error=/scicore/home/pothin/golmon00/JOB_OUT/qc_%A_%a_error.log
#SBATCH --qos=30min

module purge
module load FastQC

output_folder=$1
input_folder=$2

reads_files=(${input_folder}*.fastq.gz)
echo $reads_files

fastqc -o $1 ${reads_files[@]}
