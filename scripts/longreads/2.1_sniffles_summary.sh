#!/bin/bash -l
#SBATCH -J sniffles
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 2.1_sniffles_summary.sh DKoC53 D_RbY23_chr18_flye
# sbatch 2.1_sniffles_summary.sh DKoC53 DNeY39_chr18_flye
# sbatch 2.1_sniffles_summary.sh DLoC21 D_RbY23_chr18_flye
# sbatch 2.1_sniffles_summary.sh DLoC21 DNeY39_chr18_flye
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
cd 02_SV

for i in *_chr18_flye.vcf.gz
do
echo $i
zless $i | grep -v "#" | grep BND | wc -l 
zless $i | grep -v "#" | grep DEL | wc -l
zless $i | grep -v "#" | grep DUP | wc -l
zless $i | grep -v "#" | grep INS | wc -l
zless $i | grep -v "#" | grep INV | wc -l  
done 

# plot in long_result_results.R

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
