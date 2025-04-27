#!/bin/bash
#SBATCH --job-name="G4_lowest_4096_800"
#SBATCH --output="slurm-%j.log"
#SBATCH -p cca
##SBATCH -p preempt --qos=preempt
##SBATCH -p request --qos=request
#SBATCH -C genoa
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH -t 24:00:00    
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=you@example.com

module load modules/2.2-20230808
module load openmpi/4.0.7
module load gcc
module load hdf5/mpi-1.10.9
module load gsl/2.7
module load fftw/mpi-3.3.10

#export UCX_TLS=self,shm,ud
#export NUMBER_OF_MPI_LISTENERS_PER_NODE=2
#export UCX_LOG_LEVEL=info

cd src
make clean
make -j CONFIG=Config.sh EXEC=Gadget4
cd ..

mpiexec -n ${SLURM_NTASKS} ./Gadget4 lowest_4096_800.txt 1
