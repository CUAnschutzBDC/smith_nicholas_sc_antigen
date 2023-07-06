#! /usr/bin/env bash

#BSUB -J cellranger
#BSUB -o logs/cellranger_%J.out
#BSUB -e logs/cellranger_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

#module load cellranger/6.0.1
module load singularity/3.9.2

# Function to run snakemake
run_snakemake() {
    local num_jobs=$1
    local config_file=$2

    args='
        -o {log}.out 
        -e {log}.err 
        -J {params.job_name} 
        -R "{params.memory} span[hosts=1]"
        -n {threads} 
        -q rna '

    snakemake \
        --snakefile Snakefile \
        --drmaa "$args" \
        --jobs $num_jobs \
        --latency-wait 60 \
        --rerun-incomplete \
        --configfile $config_file \
        --singularity-args "--bind '/beevol/home/wellskri/Analysis/references/single_cell_references'" \
        --use-singularity
}

run_snakemake 12 config.yaml