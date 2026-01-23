#!/bin/bash -l
#SBATCH -J mapTE
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# cd adapterremoval/output
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/05_aDNA/1.0.1_mapTE.sh $i 00_eager_twist adapterremoval/output; done
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $2"_"$3}'); do sbatch /PATH/05_aDNA/1.0.1_mapTE.sh $i 00_eager_mybaits adapterremoval/output; done

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
#default seeding for -l is 32bp
#-n changed from 0.02 to 0.04 (stricter)

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

