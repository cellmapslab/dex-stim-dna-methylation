#!/bin/bash
#SBATCH --job-name="getRgSetFormats"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=get-rg-set-formats.out

module load R

Rscript --vanilla getRgSetFormats.R

