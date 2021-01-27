#!/bin/bash
#SBATCH --job-name="methylFiltering_quantileNMs"
#SBATCH --part=pe
#SBATCH --mem=200GB
#SBATCH --output=ms.out

module load R

Rscript --vanilla 04_getFilteredQuantNMs.R ~/github/dex-stim-dna-methylation/input_parameters.csv
