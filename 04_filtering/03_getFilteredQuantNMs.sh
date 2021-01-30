#!/bin/bash
#SBATCH --job-name="methylFiltering_quantileNMs"
#SBATCH --part=pe
#SBATCH --mem=200GB
#SBATCH --output=ms.out

module load R

SRC_DIR=~/mpip/code/dex-stim-dna-methylation

Rscript --vanilla $SRC_DIR/04_filtering/03_getFilteredQuantNMs.R $SRC_DIR/input_parameters.csv
