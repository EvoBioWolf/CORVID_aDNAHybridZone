#!/bin/bash -l
#SBATCH -J mapchr18
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 1.2_mapchr18.sh D_RbY23
# sbatch 1.2_mapchr18.sh DNeY39
# sbatch 1.2_mapchr18.sh DKoC53
# sbatch 1.2_mapchr18.sh DLoC21

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}
cd 01_mapped

samtools view -b ${1}.sorted.bam NC_046347.1 > ${1}.chr18.bam
samtools index ${1}.chr18.bam
samtools fastq -@ 8 ${1}.chr18.bam > ${1}.chr18.fastq

samtools fastq -@ 8 ${1}.unmapped.bam > ${1}.unmapped.fastq
cat ${1}.unmapped.fastq ${1}.chr18.fastq > ${1}.chr18_unmapped.fastq

pigz -p 8  ${1}.chr18_unmapped.fastq
pigz -p 8 ${1}.chr18.fastq

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
