#!/bin/bash
#SBATCH --job-name=eqtl_mapping_pl
#SBATCH --mem=200G
#SBATCH --part=pe

genotypes_trityper_dir=$1
trait_norm_filename=$2
annotation_filename=$3
coupling_filename=$4
out_mixupmapper_dir=$5

mixupmapper_dir=/home/ahryhorzhevska/mpip/tools/MixupMapper/eqtl-mapping-pipeline-1.2.4E-SNAPSHOT

java -Xmx30g -Xms30g -jar $mixupmapper_dir/eqtl-mapping-pipeline.jar --mode mixupmapper --in $genotypes_trityper_dir --out $out_mixupmapper_dir --inexp $trait_norm_filename --inexpplatform EPIC --inexpannot $annotation_filename --gte $coupling_filename
