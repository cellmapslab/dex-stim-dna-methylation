#!/bin/bash
#SBATCH --job-name="methylFiltering_quantileNBetas"
#SBATCH --part=pe
#SBATCH --mem=200GB
#SBATCH --output=betas.out

module load R

Rscript --vanilla 04_getFilteredQuantNBetas.R ~/github/dex-stim-dna-methylation/input_parameters.csv
