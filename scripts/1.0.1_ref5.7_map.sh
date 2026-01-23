#!/bin/bash -l
#SBATCH -J mapref5
#SBATCH --cpus-per-task=16
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# cd adapterremoval/output
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_ref5.7_map.sh $i 00_eager_TE adapterremoval/output; done
# cd lanemerging
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_ref5.7_map.sh $i 00_eager_WGS lanemerging; done

#sbatch 1.0.1_ref5.7_map.sh DVT014 00_eager_wgs1 lanemerging
#sbatch 1.0.1_ref5.7_map.sh DVT014 00_eager_wgs2 lanemerging
#sbatch 1.0.1_ref5.7_map.sh DVT014 00_eager_wgs3 lanemerging

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/05_aDNA"
ref="/PATH/00_pilot_2021/Corvus_cornix__S_Up_H32__v5.7_chrM.fasta"

cd $dat
cd $2/results
mkdir ext_mapped_ref5.7
cd $3

for fname in ${1}*.fq.gz
do
bwa aln -l 1024 -n 0.01 -o 2 -t 16 ${ref} ${fname} > ${dat}/${2}/results/ext_mapped_ref5.7/${1}.sai
done
#as opposed to previous seeding -l 16500, but probably does not matter as seeding is likely disabled (INT>length of sequence)
#default seeding -l 32bp

cd ${dat}/${2}/results/ext_mapped_ref5.7
for fname in ${1}.sai
do
bwa samse ${ref} ${1}.sai ${dat}/${2}/results/${3}/${1}*.fq.gz > ${1}.sam
done

#no quality threshold
samtools view -bh -@ 16 -bS ${1}.sam > ${1}_noqc.bam
samtools sort -@ 16 ${1}_noqc.bam -o ${1}_sorted_noqc.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

