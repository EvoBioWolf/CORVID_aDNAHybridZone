#!/bin/bash -l
#SBATCH -J mito
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.7_mito.sh 00_eager_TE

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
module load bcftools
module load seqtk
dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
mkdir 01_chrM
cd ${1}
cd results/deduplication

for i in *_TE
do
cd $i
# samtools index -b ${i}_rmdup.bam
# samtools view -b ${i}_rmdup.bam chrM > ${i}_chrM.bam
# samtools sort ${i}_chrM.bam > ${i}_chrM_sorted.bam
# samtools index -b ${i}_chrM_sorted.bam -@ 4
# samtools mpileup -f ${ref} -r chrM -Q 20 -q 20 -u ${i}_chrM_sorted.bam |  bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > ${dat}/01_chrM/${i}_chrM.fastq
samtools mpileup -f ${ref} -r chrM -u ${i}_chrM_sorted.bam |  bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > ${dat}/01_chrM/${i}_chrM_relaxed.fastq
# seqtk seq -aQ64 -q20 -n N ${dat}/01_chrM/${i}_chrM.fastq > ${dat}/01_chrM/${i}_chrM.fasta
seqtk seq -aQ64 -q0 -n N ${dat}/01_chrM/${i}_chrM_relaxed.fastq > ${dat}/01_chrM/${i}_chrM_relaxed.fasta
cd ../
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

