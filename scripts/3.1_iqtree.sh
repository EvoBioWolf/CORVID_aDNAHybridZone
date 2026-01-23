#!/bin/bash -l
#SBATCH -J iqtree
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 3.1_iqtree.sh enrichedallfresh_rescaled outlier_geno_maxmis_q20_plink_filtered_norm_outgroup_biallele

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}
cd $dat/02_results/iqtree
#python ${dat}/vcf2phylip/vcf2phylip.py -i ${dat}/01_angsd_${1}/${1}_${2}.vcf.gz
#iqtree -s ${1}_${2}.min4.phy -st DNA -m GTR+ASC -B 1000 -T 10
iqtree -s ${1}_${2}.min4.phy.varsites.phy -st DNA -m GTR+ASC -B 1000 -T 10

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"

