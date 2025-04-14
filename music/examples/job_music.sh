#!/bin/bash
#SBATCH --job-name="lowest"
#SBATCH --output="slurm-%j.log"
##SBATCH -p preempt --qos=preempt
##SBATCH -p cca
##SBATCH -C genoa
#SBATCH -p mem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
##SBATCH --cpus-per-task=96
#SBATCH --export=ALL
#SBATCH -t 12:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
##SBATCH --mail-user=you@email.com

module load modules/2.1.1-20230405
module load gsl
module load fftw
module load hdf5

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

./MUSIC lowest_selected_region.conf
