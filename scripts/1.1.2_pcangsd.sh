#!/bin/bash -l
#SBATCH -J pcangsd
#SBATCH --cpus-per-task=8
#SBATCH --time=10-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.2_pcangsd.sh enrichedallfresh_rescaled outlier_1111_geno_maxmis_q20_dp3
# sbatch 1.1.2_pcangsd.sh enrichedallfresh_rescaled neutral_seg72k_geno_maxmis_q20

# sbatch 1.1.2_pcangsd.sh wgs_rescaled all_geno_maxmis_q20_dp3_varonly

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"

cd ${dat}

cd ${dat}/01_angsd_${1}
pcangsd --beagle ${1}_${2}.beagle.gz --out ${1}_${2} --threads 8 \
--snp_weights -pcadapt --tree

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
