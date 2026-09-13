#!/bin/bash -l
#SBATCH -J prepgghybrid
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 4.2_prep_gghybrid.sh

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="PATH/05_aDNA"

## this script is to prep input for gghybrid ##

cd ${dat}/01_angsd_enrichedallfresh_rescaled

plink --bfile enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink --allow-extra-chr --recode structure --out ${dat}/02_results/gghybrid/enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink
plink --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink --allow-extra-chr --recode structure --out ${dat}/02_results/gghybrid/enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
