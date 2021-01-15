#!/bin/bash
#SBATCH --job-name="methylNormalization"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization.out

module load R

Rscript --vanilla 03_normalization.R ~/github/dex-stim-dna-methylation/input_parameters.csv
