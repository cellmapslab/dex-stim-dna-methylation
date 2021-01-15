#!/bin/bash
#SBATCH --job-name="methylNormalization_quantileNMs"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization_quantileNMs.out

module load R

Rscript --vanilla 03_getquantileNMs.R ~/github/dex-stim-dna-methylation/input_parameters.csv
