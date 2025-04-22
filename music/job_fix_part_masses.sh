#!/bin/bash
#SBATCH --job-name="lowest"
#SBATCH --output="/mnt/home/drennehan/ceph/trash/fix_part_masses-%j.log"
##SBATCH -p preempt --qos=preempt
#SBATCH -p cca
#SBATCH -C genoa
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --export=ALL
#SBATCH -t 3:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=you@example.com

source /mnt/home/drennehan/venvs/flathub/bin/activate

FILE_PREFIX=lowest_4096_800
START_INDEX=0
END_INDEX=7

python -u fix_part_masses.py ${FILE_PREFIX} ${START_INDEX} ${END_INDEX}
