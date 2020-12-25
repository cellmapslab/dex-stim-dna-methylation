#!/bin/bash
#SBATCH --job-name="getRgSetFormats"
#SBATCH --nodelist=hp02
#SBATCH --mem=50GB
#SBATCH --output=get-rg-set-formats.out

cd /binder/mgp/workspace/2020_DexStim_Array_Human/methylation

module load R

Rscript getRgSetFormats2.R

