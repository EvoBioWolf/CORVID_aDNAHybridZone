#!/bin/bash -l
#SBATCH -J angsd1k
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.1_angsd_1111snp.sh ancient 1111 bam_enrichedall_rescaled
# sbatch 1.1.1_angsd_1111snp.sh fresh 1111 bam_fresh


echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}/02_results
cd ./uli_pca
mkdir 01_angsd_${1}
cd 01_angsd_${1}


#hap_concensus_maxmis_q20
angsd -bam ${dat}/${3}.filelist -dohaplocall 2 -doCounts 1 -rf ../knief_snp_1based.pos \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 4 -minMinor 1 -out ${1}_${2}_hapconsensus_q20

#geno
angsd -bam ${dat}/${3}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} \
-doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -geno_minDepth 3 \
-nThreads 4 -rf ../knief_snp_1based.pos -out ${1}_${2}_geno_q20_dp3

#convert to vcf
#first created proper tfam file with proper sample name
mv ${1}_${2}_geno_q20_dp3.tfam ${1}_${2}_geno_q20_dp3_original.tfam
awk '{print $2, $1,"0","0","0","-9"}' pop.txt > ${1}_${2}_geno_q20_dp3.tfam
plink --tfile ${1}_1111_geno_q20_dp3 --allow-no-sex --allow-extra-chr -make-bed --out ${1}_1111_geno_q20_dp3_plink
plink --bfile ${1}_1111_geno_q20_dp3_plink --allow-extra-chr --recode vcf --out ${1}_1111_geno_q20_dp3_plink
plink --bfile ${1}_1111_geno_q20_dp3_plink --allow-extra-chr --pca --out ${1}_1111_geno_q20_dp3_plink

#MODIFY OUTPUT
zcat ${1}_1111_geno_q20_dp3.geno.gz | \
awk '{printf "%s", $1":"$2; for (i=3; i<=NF; i++) {printf "\t%s/%s", substr($i, 1, 1), substr($i, 2, 1);}print "";}' | \
datamash transpose > ${1}_1111_geno_q20_dp3.geno.transposed.txt

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
