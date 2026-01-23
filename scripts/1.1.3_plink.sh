#!/bin/bash -l
#SBATCH -J plink
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled outlier_geno_maxmis_q20
# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_trans
# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled outlier_1111_geno_maxmis_q20

# sbatch 1.1.3_plink.sh wgs_rescaled all_geno_maxmis_q20_dp3_varonly
# sbatch 1.1.3_plink.sh wgs_rescaled all_geno_maxmis_q20_dp5_varonly

# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled neutral_seg72k_geno_maxmis_q20
# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_dp3
# sbatch 1.1.3_plink.sh enrichedallfresh_rescaled outlier_1111_geno_maxmis_q20_dp3

# sbatch 1.1.3_plink.sh wgs_rescaled_outgroup all_geno_maxmis_q20_dp3

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
#using version 1.9

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd 01_angsd_${1}

#create proper tfam file
mv ${1}_${2}.tfam ${1}_${2}_original.tfam
awk '{print $2, $1, "0", "0", "0", "-9"}' pop.txt >> ${1}_${2}.tfam

plink --tfile ${1}_${2} --allow-no-sex --allow-extra-chr -make-bed --out ${1}_${2}_plink

plink --bfile ${1}_${2}_plink --allow-extra-chr --recode vcf --out ${1}_${2}_plink

plink --bfile ${1}_${2}_plink --allow-extra-chr --pca --out ${1}_${2}_plink

plink --bfile ${1}_${2}_plink --allow-extra-chr --missing --out ${1}_${2}_plink

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
