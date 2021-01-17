#!/bin/bash
#SBATCH --job-name="methylNormalization_quantileNBetas"
#SBATCH --part=pe
#SBATCH --mem=300GB
#SBATCH --output=methyl_normalization_quantileNBetas.out

module load R

Rscript --vanilla 03_getquantileNBetas.R ~/github/dex-stim-dna-methylation/input_parameters.csv
