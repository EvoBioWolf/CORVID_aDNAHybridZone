#!/bin/bash -l
#SBATCH -J align
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 3.2_align.sh D_RbY23
# sbatch 3.2_align.sh DNeY39
# sbatch 3.2_align.sh DKoC53
# sbatch --dependency=afterok:3267187 3.2_align.sh DKoC45
# sbatch 3.2_align.sh DKoC21

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}

#awk '/^>NC_046347.1/{flag=1} /^>/{if(flag && !/^>NC_046347.1/){flag=0}} flag' GCF_000738735.6_ASM73873v6_genomic.fna > ASM73873v6_chr18.fasta
#minimap2 -d ${dat}/ASM73873v6_chr18.mmi ${dat}/ASM73873v6_chr18.fasta

cd 03_assemblies

# See how Flye contigs map along the reference
# minimap2 -ax asm5 ${dat}/ASM73873v6_chr18.mmi ${dat}/03_flye/${1}/${1}_chr18_flye.fasta | samtools sort -o ${1}_chr18_flye_ref.bam
# samtools index ${1}_chr18_flye_ref.bam
# qualimap bamqc -bam ${1}*_flye_ref.bam -c --java-mem-size=32G


# Piece together tha contigs according to reference but risk structural bias
# ragtag.py scaffold -o ${1}_chr18_flye_ragtag ${dat}/ASM73873v6_chr18.fasta ${dat}/03_flye/${1}/${1}_chr18_flye.fasta

cd ${dat}/03_assemblies/${1}_chr18_flye_ragtag
# awk '/^>NC_046347.1/{f=1; print; next} /^>/{f=0} f' ragtag.scaffold.fasta > ${1}_chr18_flye_ragtag.fasta
# minimap2 -d ${1}_chr18_flye_ragtag.mmi ${1}_chr18_flye_ragtag.fasta
# minimap2 -ax asm5 --eqx ${dat}/ASM73873v6_chr18.mmi ${1}_chr18_flye_ragtag.fasta | samtools sort -o ${1}_chr18_flye_ragtag_ref.bam
# samtools index ${1}_chr18_flye_ragtag_ref.bam
# minimap2 -ax asm5 --eqx ${1}_chr18_flye_ragtag.mmi ${dat}/ASM73873v6_chr18.fasta | samtools sort -o ref_to_${1}_chr18_flye_ragtag.bam
# samtools index ref_to_${1}_chr18_flye_ragtag.bam
# qualimap bamqc -bam ref_to_${1}_chr18_flye_ragtag.bam -c --java-mem-size=32G

#creating a mapping loop
# minimap2 -ax asm5 --eqx ../D_RbY23_chr18_flye_ragtag/D_RbY23_chr18_flye_ragtag.mmi ${1}_chr18_flye_ragtag.fasta | samtools sort -o ${1}_to_D_RbY23_chr18_flye_ragtag.bam
# samtools index ${1}_to_D_RbY23_chr18_flye_ragtag.bam
# minimap2 -ax asm5 --eqx ../DNeY39_chr18_flye_ragtag/DNeY39_chr18_flye_ragtag.mmi ${1}_chr18_flye_ragtag.fasta | samtools sort -o ${1}_to_DNeY39_chr18_flye_ragtag.bam
# samtools index ${1}_to_DNeY39_chr18_flye_ragtag.bam
# minimap2 -ax asm5 --eqx ../DKoC53_chr18_flye_ragtag/DKoC53_chr18_flye_ragtag.mmi ${1}_chr18_flye_ragtag.fasta | samtools sort -o ${1}_to_DKoC53_chr18_flye_ragtag.bam
# samtools index ${1}_to_DKoC53_chr18_flye_ragtag.bam
minimap2 -ax asm5 --eqx ../DKoC21_chr18_flye_ragtag/DKoC21_chr18_flye_ragtag.mmi ${1}_chr18_flye_ragtag.fasta | samtools sort -o ${1}_to_DKoC21_chr18_flye_ragtag.bam
samtools index ${1}_to_DKoC21_chr18_flye_ragtag.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
