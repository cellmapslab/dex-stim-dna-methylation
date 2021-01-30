#!/bin/bash
#SBATCH --job-name="methylNormalization"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization.out

module load R

SRC_DIR=~/mpip/code/dex-stim-dna-methylation

Rscript --vanilla $SRC_DIR/03_normalization/01_normalization.R $SRC_DIR/input_parameters.csv
