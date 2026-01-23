#!/bin/bash -l
#SBATCH -J regele
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 3.3_regelement.sh

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}
cd $dat/02_results/uli_pca
cd regulatory_element

# I created a files to look for regulatory elements 
# These are 100bp regions upstream and downstream of the SNPs potentially in regulatory regions (SNPS just up/ down of a gene)
# So each target region is 212bp, including the the SNP in the middle (bed:start=-101fromSNPpos, end=+100fromSNPpos)
# Since gff is restricted to intergenic/ upstream/ downstream annotation - not informative enough
# Download known vertebrate motifs from JASPAR then run Fimo on MEME

bedtools getfasta -fi ${ref} -bed region_0based.bed -s -name -fo upstream_downstream_regions.fa
#fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt upstream_downstream_regions.fa

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"

