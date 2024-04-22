#! /usr/bin/env bash

#BSUB -J immcantation
#BSUB -o logs/immcantation_%J.out
#BSUB -e logs/immcantation_%J.err
#BSUB -R "select[mem>100] rusage[mem=100]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

#module load cellranger/6.0.1
module load singularity/3.9.2

base_dir=/beevol/home/wellskri/Analysis

smith_dir=${base_dir}/Mia_Smith/Catherine_Nicolas/full_antigen_pos_data
nakayama_dir=${base_dir}/Maki_Nakayama/pln_scrna_seq

singularity_container=${smith_dir}/docker/immcantation/immcantation.sif

sample1=${smith_dir}/results/R_analysis/merged/immcantation_combined.tsv
sample2=${nakayama_dir}/results/R_analysis/merged/immcantation_combined.tsv

combined=${smith_dir}/results/R_analysis/merged/immcantation_combined_pln.tsv

all_output=$smith_dir/results/R_analysis/merged/define_clones/immcantation_combined_pln-pass.tsv

awk 'NR==1 {print; next} FNR>1' $sample1 $sample2 > $combined

# Run immcantation define clones
singularity exec DefineClones.py \
	-d $combined \
	-o $all_output \
	--nproc 10 \
	--model ham \
	--norm len \
	--dist 0.15