#!/bin/bash

#SBATCH --job-name=benchmarkDLPFC
#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --time=60:00:00
#SBATCH --mem=24000
#SBATCH --mail-user=ezhao@fredhutch.org
#SBATCH --mail-type=END,FAIL
ml R/3.6.2-foss-2016b-fh1

cd /fh/fast/gottardo_r/ezhao_working/DLPFC/
R CMD BATCH --no-save --no-restore ./benchmarkDLPFC.R benchmarkDLPFC_${SLURM_ARRAY_TASK_ID}.Rout
