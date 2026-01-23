#!/bin/bash -l
#SBATCH -J mapWGS
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# cd lanemerging
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh $i 00_eager_wgs1 lanemerging; done
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh $i 00_eager_wgs2 lanemerging; done
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh $i 00_eager_wgs3 lanemerging; done

# sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh KCZ003_UDG 00_eager_wgs1 adapterremoval/output
# sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh NCP001_UDG 00_eager_wgs2 adapterremoval/output
# sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh WMP006_WGS 00_eager_wgs2 adapterremoval/output
# sbatch /PATH/05_aDNA/1.0.1_mapWGS.sh NCP001_UDG 00_eager_wgs3 adapterremoval/output

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd ${2}/results
mkdir ext_mapped
cd ${3}

echo ${1}

for fname in *${1}*.fq.gz
do
bwa aln -l 1024 -n 0.04 -o 2 -t 8 ${ref} ${fname} > ${dat}/${2}/results/ext_mapped/${1}.sai
done

cd ${dat}/${2}/results/ext_mapped
for fname in ${1}.sai
do
bwa samse ${ref} ${1}.sai ${dat}/${2}/results/${3}/*${1}*.fq.gz > ${1}.sam
done

#no quality threshold
samtools view -bh -@ 8 -bS ${1}.sam > ${1}_noqc.bam
samtools sort -@ 8 ${1}_noqc.bam -o ${1}_sorted_noqc.bam
rm ${1}_noqc.bam

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

