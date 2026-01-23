#!/bin/bash -l
#SBATCH -J mitoalign
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.7_mito_align.sh 1000

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
module load mafft
module load raxml

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd 01_chrM

# for fname in *_chrM.fasta
# do
# base=${fname%_chrM.fasta*}
# grep -o N ${fname} | wc -l >> wc.tmp 
# done

# ls *_chrM.fasta > wc2.tmp
# paste wc2.tmp wc.tmp > 01_missing_count.txt
# rm wc*.tmp

# for fname in *_chrM_relaxed.fasta
# do
# base=${fname%_chrM_relaxed.fasta*}
# grep -o N ${fname} | wc -l >> wc_relaxed.tmp
# done

# ls *_chrM_relaxed.fasta > wc2_relaxed.tmp
# paste wc2_relaxed.tmp wc_relaxed.tmp > 01_missing_count_relaxed.txt
# rm wc*_relaxed.tmp

#select sequences with less than 8000 Ns (<50% missingness)
# for fname in $(awk '{if ($2<8000) print $1}' 01_missing_count.txt)
# do
# base=${fname%_TE_chrM.fasta*}
# awk 'NR >1' ${base}_TE_chrM.fasta | awk -v a=${base} 'BEGIN {print ">"a} {print $0}' >> allmito_selected.fasta
# done
# mafft --localpair --thread 8 --maxiterate 1000 allmito_selected.fasta > allmito_selected_linsi_aligned.fasta

#add correctly aligned genbank corvids
# cat genbank_realigned.fasta H32_chrM.fasta allmito_selected.fasta > allmito_selected_outgroup.fasta
# mafft --localpair --thread 8 --maxiterate 1000 allmito_selected_outgroup.fasta > allmito_selected_outgroup_linsi_aligned.fasta

# raxmlHPC -T 8 -s allmito_selected_outgroup_linsi_aligned.fasta -m GTRGAMMA -N ${1} -n allmito_selected_outgroup_linsi_aligned_bs${1} -p 12367 -f a -x 16789

#all samples
for fname in $(awk '{if ($2<16000) print $1}' 01_missing_count.txt)
do
base=${fname%_TE_chrM.fasta*}
awk 'NR >1' ${base}_TE_chrM.fasta | awk -v a=${base} 'BEGIN {print ">"a} {print $0}' >> allmito.fasta
done
mafft --localpair --thread 8 --maxiterate 1000 allmito.fasta > allmito_linsi_aligned.fasta
cat genbank_realigned.fasta H32_chrM.fasta allmito.fasta > allmito_outgroup.fasta
mafft --localpair --thread 8 --maxiterate 1000 allmito_outgroup.fasta > allmito_outgroup_linsi_aligned.fasta
raxmlHPC -T 8 -s allmito_outgroup_linsi_aligned.fasta -m GTRGAMMA -N ${1} -n allmito_outgroup_linsi_aligned_bs${1} -p 12367 -f a -x 16789

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

