#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-3:00:00
#SBATCH --mem=24G

WD=/cfs/klemming/scratch/e/ebazzica/selscans
THR=4
BAYP=${WD}/baypass_public-v3.1/sources/g_baypass
ACFILES=${WD}/subsets_acs
OUTDIR=${WD}/baypass_core

subset=${1}
acfile=${ACFILES}/60.Neff.${subset}.ac
outpref=${OUTDIR}/${subset}

${BAYP} -npop 53 -gfile ${acfile} -outprefix ${outpref} -nthreads ${THR}
