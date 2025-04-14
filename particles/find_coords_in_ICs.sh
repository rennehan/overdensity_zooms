#!/bin/bash -l
#########################################################
#SBATCH -J find_coords_in_IC
#SBATCH -p cca
##SBATCH -p preempt --qos=preempt
#SBATCH -C rome
##SBATCH --mail-user=you@email.com
#SBATCH --mail-type=FAIL 
#SBATCH -o ./find_coords_in_IC-%j.out
#########################################################
#SBATCH --time=6:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --array=0-63
#########################################################

source /mnt/home/drennehan/venvs/flathub/bin/activate

DATA_DIR=/mnt/home/drennehan/ceph/simulating_the_universe/overdensity_intersections/sims/4096_800/output
NUM_FILES=1
SNAPSHOT_IDX=111
ICS_FILE=/mnt/home/drennehan/ceph/simulating_the_universe/overdensity_intersections/boxes/4096_800/4096_800.${SLURM_ARRAY_TASK_ID}.hdf5
INDEX_LIST_FILE=${DATA_DIR}/lowest_cell_pids_111.pkl
BOX_SIZE_IN_KPC=800000.0
IC_FILE_IDX=${SLURM_ARRAY_TASK_ID}

python -u find_coords_in_ICs.py ${DATA_DIR} ${NUM_FILES} ${SNAPSHOT_IDX} ${ICS_FILE} ${INDEX_LIST_FILE} ${BOX_SIZE_IN_KPC} ${IC_FILE_IDX}

