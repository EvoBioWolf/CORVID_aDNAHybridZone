#!/bin/bash -l
#SBATCH -J remap
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 1.3_remap.sh DNeY39 DKoC53_chr18_flye_ragtag
# sbatch 1.3_remap.sh D_RbY23 DKoC53_chr18_flye_ragtag
# sbatch 1.3_remap.sh DKoC53 D_RbY23_chr18_flye
# sbatch 1.3_remap.sh DKoC53 DNeY39_chr18_flye
# sbatch 1.3_remap.sh DLoC21 D_RbY23_chr18_flye
# sbatch 1.3_remap.sh DLoC21 DNeY39_chr18_flye
# sbatch 1.3_remap.sh DNeY39 DLoC21_chr18_flye
# sbatch 1.3_remap.sh D_RbY23 DLoC21_chr18_flye
# sbatch 1.3_remap.sh DNeY39 DKoC53_chr18_flye
# sbatch 1.3_remap.sh D_RbY23 DKoC53_chr18_flye


echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"
#newrf="/PATH/06_longreads/03_assemblies/DKoC53_chr18_flye_ragtag/DKoC53_chr18_flye_ragtag.mmi"
#newrf="/PATH/06_longreads/03_assemblies/DNeY39_chr18_flye_ragtag/DNeY39_chr18_flye_ragtag.mmi"
#newrf="/PATH/06_longreads/03_flye/D_RbY23/D_RbY23_chr18_flye"
#newrf="/PATH/06_longreads/03_flye/DNeY39/DNeY39_chr18_flye"
#newrf="/PATH/06_longreads/03_flye/DKoC21/DKoC21_chr18_flye"
newrf="/PATH/06_longreads/03_flye/DKoC53/DKoC53_chr18_flye"

cd ${dat}

cd 02_remap
minimap2 -ax map-ont -t 8 ${newrf}.fasta ${dat}/01_mapped/${1}.chr18.fastq.gz | samtools view -b -@ 8 | samtools sort -@ 8 -o ${1}_ref${2}.remap.bam
samtools index ${1}_ref${2}.remap.bam

#cd 01_mapped
#minimap2 -d ${newrf}.mmi ${newrf}.fasta
#minimap2 -ax map-ont -t 8 ${newrf}.mmi ${1}.chr18_unmapped.fastq.gz | samtools view -b -@ 8 | samtools sort -@ 8 -o ../02_remap/${1}_ref${2}.remap.bam
#samtools index ../02_remap/${1}_ref${2}.remap.bam

#qualimap bamqc -bam ${1}.remap.bam -c --java-mem-size=32G

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
