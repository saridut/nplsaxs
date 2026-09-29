#!/usr/bin/env bash

#SBATCH --job-name=saxsds # Jobname
#SBATCH --partition=Cascade
#SBATCH --ntasks=1      
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
##SBATCH --hint=nomultithread
#SBATCH --array=0-5 ##Should go to 500
#SBATCH --output=%x.o%j      
#SBATCH --error=%x.o%j      
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=benjamin.abecassis@ens-lyon.fr
#SBATCH --time=1-00:00:00

#module purge
#module use /applis/PSMN/debian13/Cascade/modules/all
#module load foss/2025a
#conda activate saxs

SLURM_ARRAY_TASK_ID=4

dn='outdir'
if [ ! -d "${dn}" ]; then
    mkdir -p ${dn}
fi

length=100  #nm
width=10    #nm
ML=3        #num monolayers

#python -m nplsaxs.cli create \
python -m nplsaxs.cli create \
    --rdist lognormal \
    --rmin 5 \
    --rmax 20 \
    --rpd 1.2 \
    --pdist lognormal \
    --pmin 10 \
    --pmax 50 \
    --ppd 1.2 \
    --outfile "${dn}/ds_${SLURM_ARRAY_TASK_ID}.npz" \
    --nsamp 20 \
    --phi 0.001 \
    --npart 512 \
    --calculator 'AESDebye' \
    --pattern_type 'x'      \
    --Qbeg 0.01 \
    --Qend 80    \
    --Qstep 0.02 \
    --nthreads -1 \
     ${length} ${width} ${ML}
