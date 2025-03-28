#!/bin/bash
#
#SBATCH --array=0-2
#SBATCH --cpus-per-task=48
#SBATCH --job-name=neotest
#SBATCH --output=slurm_%a.out
#SBATCH --mem=96000
#SBATCH --time=50:00:00
#SBATCH --partition=ceewater_cjgleason-cpu
#SBATCH --error=slurm.%N.%j.err
/work/pi_cjgleason_umass_edu/.conda/envs/lightweight/lib/R/bin/Rscript --vanilla slurm_run.R
