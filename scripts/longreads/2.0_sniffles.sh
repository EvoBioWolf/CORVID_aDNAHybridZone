#!/bin/bash -l
#SBATCH -J sniffles
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 2.0_sniffles.sh D_RbY23

# sbatch 2.0_sniffles.sh D_RbY23 DKoC53
# sbatch 2.0_sniffles.sh DNeY39 DKoC53
# sbatch 2.0_sniffles.sh DKoC53 DNeY39
# sbatch 2.0_sniffles.sh DKoC53 D_RbY23_chr18_flye
# sbatch 2.0_sniffles.sh DKoC53 DNeY39_chr18_flye
# sbatch 2.0_sniffles.sh DLoC21 D_RbY23_chr18_flye
# sbatch 2.0_sniffles.sh DLoC21 DNeY39_chr18_flye
# sbatch 2.0_sniffles.sh D_RbY23 DKoC53_chr18_flye
# sbatch 2.0_sniffles.sh DNeY39 DKoC53_chr18_flye
# sbatch 2.0_sniffles.sh D_RbY23 DLoC21_chr18_flye
# sbatch 2.0_sniffles.sh DNeY39 DLoC21_chr18_flye

echo $(date)
STARTTIME=$(date +%s)

conda activate nanoplot #changed, nanoplot removed
dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}

cd 02_remap
sniffles --input ${1}_ref${2}.remap.bam --vcf ${dat}/02_SV/${1}_ref${2}.vcf.gz --snf ${dat}/02_SV/${1}_ref${2}.snf -t 8

sniffles --input ${1}_ref${2}.remap.bam --vcf ${dat}/02_SV/${1}_ref${2}_min5.vcf.gz --minsupport 5 --snf ${dat}/02_SV/${1}_ref${2}_min5.snf -t 8 

#cd 01_mapped
#sniffles --input ${1}.chr18.bam --vcf ${dat}/02_SV/${1}.vcf.gz --snf ${dat}/02_SV/${1}.snf -t 8

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
