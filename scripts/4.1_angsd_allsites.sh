#!/bin/bash -l
#SBATCH -J angsd18
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 4.1_angsd_allsites.sh enriched_allsites
# sbatch 4.1_angsd_allsites.sh enriched_allsites_carrion
# sbatch 4.1_angsd_allsites.sh enriched_allsites_hooded
echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
#ref="/PATH/06_longreads/ASM73873v6_chr18.fasta"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic_renamed.fna"

cd ${dat}/01_angsd_${1}

#hap_concensus_maxmis_q20 (not ran)
# angsd -bam bam_${1}.filelist -dohaplocall 2 -doCounts 1 \
# -checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
# -nThreads 8 -minMinor 2 -setMinDepthInd 3 -out ${1}_hapconsensus_q20_dp3
# then run haploidconversion.ipynb to convert to vcf 

#geno no dp cut-off (not ran)
# angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
# -checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -SNP_pval 1e-6 \
# -nThreads 8 -out ${1}_geno_q20

# mv ${1}_geno_q20.tfam ${1}_geno_q20_original.tfam
# awk '{print $2, $1, "0", "0", "0", "-9"}' pop.txt >> ${1}_geno_q20.tfam
# plink --tfile ${1}_geno_q20 --allow-no-sex --allow-extra-chr -make-bed --out ${1}_geno_q20
# plink --bfile ${1}_geno_q20 --allow-extra-chr --recode vcf --out ${1}_geno_q20
# plink --bfile ${1}_geno_q20 --allow-extra-chr --pca --out ${1}_geno_q20
# plink --bfile ${1}_geno_q20 --allow-extra-chr --missing --out ${1}_geno_q20
# bgzip ${1}_geno_q20.vcf
# tabix -p vcf ${1}_geno_q20.vcf.gz

#geno present
# angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
# -checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -SNP_pval 1e-6 \
# -nThreads 16 -geno_minDepth 3 -out ${1}_geno_q20_dp3

# mv ${1}_geno_q20_dp3.tfam ${1}_geno_q20_dp3_original.tfam
awk '{print $2, $1, "0", "0", "0", "-9"}' pop.txt >> ${1}_geno_q20_dp3.tfam
plink --tfile ${1}_geno_q20_dp3 --allow-extra-chr --chr-set 28 -make-bed --out ${1}_geno_q20_dp3
plink --bfile ${1}_geno_q20_dp3 --allow-extra-chr --chr-set 28 --recode vcf --out ${1}_geno_q20_dp3
plink --bfile ${1}_geno_q20_dp3 --allow-extra-chr --chr-set 28 --pca --out ${1}_geno_q20_dp3
plink --bfile ${1}_geno_q20_dp3 --allow-extra-chr --chr-set 28 --missing --out ${1}_geno_q20_dp3
sed 's/^\(##contig=<ID=\)chr/\1/' ${1}_geno_q20_dp3.vcf | sed 's/^chr//' | bgzip > ${1}_geno_q20_dp3.vcf.gz
tabix -p vcf ${1}_geno_q20_dp3.vcf.gz
rm ${1}_geno_q20_dp3.vcf

# split vcf into per chr
for CHR in {1..28} 1A 4A
do
bcftools view -r ${CHR} -Oz -o ${1}_geno_q20_dp3_chr${CHR}.vcf.gz ${1}_geno_q20_dp3.vcf.gz
bcftools index ${1}_geno_q20_dp3_chr${CHR}.vcf.gz
done

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
