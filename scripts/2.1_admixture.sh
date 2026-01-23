#!/bin/bash -l
#SBATCH -J admix
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# for i in {1..10}; do sbatch 2.1_admixture.sh enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink $i poplist 01_angsd_enrichedallfresh_rescaled; done
# for i in {1..10}; do sbatch 2.1_admixture.sh enrichedallfresh_rescaled_neutral_geno_maxmis_q20_trans_plink $i poplist 01_angsd_enrichedallfresh_rescaled; done

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

echo $(date)
STARTTIME=$(date +%s)

cd $dat/$4
#sed -i "s/scaffold_//g" ${1}.bim
admixture --cv ${1}.bed $2 -j4 | tee log${2}.out 

cd ${dat}/02_results/admixture
mkdir $1
cd ./$1
mv $dat/$4/log${2}.out ./
mv $dat/$4/${1}.${2}.Q ./
mv $dat/$4/${1}.${2}.P ./

conda activate renv
Rscript ../admixtureplot.R $1 $2 ${dat}/02_results/admixture/$1 ../$3.txt

#when all Ks have completed
# grep -h CV log*.out
# conda activate base
# convert -append *.pdf all_merged.pdf

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this 
task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this 
task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete 
this task"
