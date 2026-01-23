#!/bin/bash -l
#SBATCH -J flye
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 3.0_flye.sh D_RbY23
# sbatch 3.0_flye.sh DNeY39
# sbatch 3.0_flye.sh DKoC53
# sbatch 3.0_flye.sh DLoC21


echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}
cd 03_flye

mkdir ${1}
flye --nano-hq ${dat}/01_mapped/${1}.chr18.fastq.gz \
-g 12m -o ${1} \
-t 8 --min-overlap 5000
cd ${1}
rename assembly ${1}_chr18_flye_2 *

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
