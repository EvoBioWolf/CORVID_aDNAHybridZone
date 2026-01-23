#!/bin/bash -l
#SBATCH -J emu
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.3_emu.sh

echo $(date)
STARTTIME=$(date +%s)

conda activate emu

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}

#install emu
#git clone https://github.com/Rosemeis/emu.git
#conda env create -f emu/environment.yml

cd 01_angsd_enrichedallfresh_rescaled

#neutral
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_trans_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_trans_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink.emu

#outlier
emu --bfile enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_plink.emu
emu --bfile enrichedallfresh_rescaled_outlier_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_outlier_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_geno_maxmis_q20_plink.emu


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
