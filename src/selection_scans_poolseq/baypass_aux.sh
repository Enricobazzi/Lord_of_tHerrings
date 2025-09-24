#!/bin/bash -l
#SBATCH -A naiss2025-22-1002
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 0-4:00:00
#SBATCH --mem=24G

WD=/cfs/klemming/scratch/e/ebazzica/selscans
THR=4
BAYP=${WD}/baypass_public-v3.1/sources/g_baypass
ACFILES=${WD}/subsets_acs
OMEGAS=${WD}/baypass_core
EFILE=${WD}/env_files/sst_mean.txt
OUTDIR=${WD}/baypass_aux

subset=${1}
acfile=${ACFILES}/60.Neff.${subset}.ac
omega=${OMEGAS}/${subset}_mat_omega.out
outpref=${OUTDIR}/sst_${subset}

${BAYP} -npop 53 -gfile ${acfile} -efile ${EFILE} -auxmodel -scalecov -omegafile ${omega} -outprefix ${outpref} -nthreads ${THR}
