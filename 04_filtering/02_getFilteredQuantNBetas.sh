#!/bin/bash
#SBATCH --job-name="methylFiltering_quantileNBetas"
#SBATCH --part=pe
#SBATCH --mem=200GB
#SBATCH --output=betas.out

module load R

SRC_DIR=~/mpip/code/dex-stim-dna-methylation

Rscript --vanilla $SRC_DIR/04_filtering/02_getFilteredQuantNBetas.R $SRC_DIR/input_parameters.csv
