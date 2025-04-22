#!/bin/bash -l
#########################################################
#SBATCH -J generate_overdensity_cells
#SBATCH -p cca
##SBATCH -p preempt --qos=preempt
##SBATCH -p mem
##SBATCH --mail-user=you@email.com
#SBATCH --mail-type=FAIL 
#SBATCH -o /mnt/home/drennehan/ceph/trash/generate_overdensity_cells-%j.out
#########################################################
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --array=26,43,58,72,132,133,135,138,141,144,147,151,156,163
#SBATCH --array=44-57,59-71,73-110
#SBATCH --cpus-per-task=96
#########################################################
# z     snap_idx
# 6     111
# 9     72
# 11    58
# 14    43
# 20    26

source /mnt/home/drennehan/venvs/flathub/bin/activate

DATA_DIR=/mnt/home/drennehan/ceph/simulating_the_universe/overdensity_intersections/sims/4096_800/output
NUM_FILES=8
SNAPSHOT_IDX=${SLURM_ARRAY_TASK_ID}
FILE_TYPE="gadget"      # gadget or swift

python -u generate_overdensity_cells.py ${DATA_DIR} ${NUM_FILES} ${SNAPSHOT_IDX} ${FILE_TYPE}

#########################################################
