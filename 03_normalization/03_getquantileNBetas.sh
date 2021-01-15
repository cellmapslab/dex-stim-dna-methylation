#!/bin/bash
#SBATCH --job-name="methylNormalization-quantileNMs"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization-quantileNMs.out

module load R

Rscript --vanilla 03_getquantileNBetas.R ~/github/dex-stim-dna-methylation/input_parameters.csv
