#!/bin/bash
#SBATCH --job-name=cali_plot
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --output=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/job_out/ampSeq.LOG
#SBATCH --error=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/job_out/ampSeq.o
#SBATCH --qos=1day

########################
# Submit the ampSeq analysis to the scicore cluster
# 109.12.2022
# monica.golumbeanu@unibas.ch
#
# To use this script, run the command on the scicore terminal:
# sbatch submit_to_cluster.sh input_folder output_folder
# where input_folder needs to contain all the files needed for the analysis
# (a description of this files is provided in the wiki:
# https://git.scicore.unibas.ch/swisstph-malaria-genotyping/ampseq_pipeline/-/wikis/home)
# output_folder is the folder where the output will be saved.
#######################

input_folder=$1
output_folder=$2

module load R/4.1.0-foss-2018b

echo "Running AmpSeq analysis ..."

Rscript haplotypR_analysis.R $input_folder $output_folder

echo "Finished."

