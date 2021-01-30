#!/bin/bash
#SBATCH --job-name="methylNormalization_quantileNMs"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization_quantileNMs.out

module load R

SRC_DIR=~/mpip/code/dex-stim-dna-methylation

Rscript --vanilla $SRC_DIR/03_normalization/03_getquantileNMs.R $SRC_DIR/input_parameters.csv
