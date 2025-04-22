#!/bin/bash -l
#########################################################
#SBATCH -J find_pids_in_cell
##SBATCH -p cca
##SBATCH -p preempt --qos=preempt
#SBATCH -p mem
##SBATCH --mail-user=you@email.com
#SBATCH --mail-type=FAIL 
#SBATCH -o /mnt/home/drennehan/ceph/trash/find_pids_in_cell-%j.out
#########################################################
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=26,43,58,72
#SBATCH --cpus-per-task=96
#########################################################

source /mnt/home/drennehan/venvs/flathub/bin/activate

DATA_DIR=/mnt/home/drennehan/ceph/simulating_the_universe/overdensity_intersections/sims/4096_800/output
NUM_FILES=8
SNAPSHOT_IDX=${SLURM_ARRAY_TASK_ID}
FILE_TYPE="gadget"      # gadget or swift

python -u find_particles_in_cell.py ${DATA_DIR} ${NUM_FILES} ${SNAPSHOT_IDX} ${FILE_TYPE}

#########################################################
