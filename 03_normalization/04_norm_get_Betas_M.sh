#!/bin/bash
#SBATCH --job-name="methylNormalization-get-betas-m"
#SBATCH --part=pe
#SBATCH --mem=50GB
#SBATCH --output=methyl_normalization-qet-betas-m.out

sbatch --dependency=afterany:2589749 03_getquantileNBetas.sh
sbatch --dependency=afterany:2589749 03_getquantileNMs.sh
