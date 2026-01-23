#!/bin/bash -l
#SBATCH -J mapref5fresh
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

#sbatch 1.0.1_ref5.7_fresh.sh

conda activate biotools 
#bwa=0.7.17
#samtools=1.7-1
#picard.jar=2.25.7

ref="/PATH/00_pilot_2021/Corvus_cornix__S_Up_H32__v5.7_chrM.fasta"
dat="/PATH/05_aDNA"

echo $(date)
STARTTIME=$(date +%s)
cd $dat

# adapter trim
#trim adapter using cutadapt

cd fresh_rawreads
for fname in *pair_1*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R1.fq.gz
done

for fname in *pair_2*.fq.gz
do
base=${fname%.fq.gz*}
mv ${base}.fq.gz ${base}.R2.fq.gz
done

rename "pair_1_Illumina" "Illumina" *
rename "pair_2_Illumina" "Illumina" *

for fname in *.R1.fq.gz
do
base=${fname%.R1.fq.gz*}
cutadapt -j 8 -m 30 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ../fresh_trimmed/${base}_trimmed.R1.fq.gz -p ../fresh_trimmed/${base}_trimmed.R2.fq.gz ${base}.R1.fq.gz ${base}.R2.fq.gz
done

echo bwamem
cd ../fresh_trimmed
for fname in *_trimmed.R1.fq.gz
do
base=${fname%_trimmed.R1.fq.gz*}
bwa mem -t 16 -M ${ref} ${base}_trimmed.R1.fq.gz ${base}_trimmed.R2.fq.gz > ../fresh_mapped/${base}.sam 
samtools view -bhq 6 -@ 16 ../fresh_mapped/${base}.sam > ../fresh_mapped/${base}.bam
samtools sort -m 4G -@ 16 ../fresh_mapped/${base}.bam -o ../fresh_mapped/${base}.sorted.bam
done

cd ../fresh_mapped
mkdir bam_sorted sam bam
mv *.sorted.bam ./bam_sorted
mv *.sam ./sam
mv *.bam ./bam

#rm *.sam
#rm *.bam
# then with continue pipeline with eager2

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) sec to complete this task"
HOUR=3600
echo "It takes $((($ENDTIME - $STARTTIME) / $HOUR)) hr to complete this task"

