#!/bin/bash
#SBATCH --job-name="methylNormalization_quantileNBetas"
#SBATCH --part=pe
#SBATCH --mem=300GB
#SBATCH --output=methyl_normalization_quantileNBetas.out

module load R

SRC_DIR=~/mpip/code/dex-stim-dna-methylation

Rscript --vanilla 02_getquantileNBetas.R $SRC_DIR/input_parameters.csv
