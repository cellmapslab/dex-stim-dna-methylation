#!/bin/bash
#SBATCH --job-name="methylQC"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_qc.out

module load R

Rscript --vanilla 01_qc.R
