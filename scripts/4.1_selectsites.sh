#!/bin/bash -l
#SBATCH -J select18
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 4.1_selectsites.sh enriched_chr18 

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/06_longreads/ASM73873v6_chr18.fasta"

cd ${dat}/01_angsd_${1}

#44samples in vcf, 0.03: 1 missing sample allowed
#bcftools view -i 'F_MISSING=0' ${1}_geno_q20_dp3.vcf.gz -Oz -o ${1}_geno_q20_dp3_miss0.vcf.gz
bcftools view -i 'F_MISSING=0.03' ${1}_geno_q20_dp3.vcf.gz -Oz -o ${1}_geno_q20_dp3_miss0.03.vcf.gz

bcftools view -i 'F_MISSING=0' ${1}_geno_q20.vcf.gz -Oz -o ${1}_geno_q20_miss0.vcf.gz
bcftools view -i 'F_MISSING=0.03' ${1}_geno_q20.vcf.gz -Oz -o ${1}_geno_q20_miss0.03.vcf.gz


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
