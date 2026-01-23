#!/bin/bash -l
#SBATCH -J variant
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.3_wgs_variant.sh wgs_rescaled all_geno_maxmis_q20_plink 
# sbatch 1.3_wgs_variant.sh wgs_rescaled all_geno_maxmis_q20_dp3_varonly_plink
# sbatch 1.3_wgs_variant.sh wgs_rescaled_outgroup all_geno_maxmis_q20_dp3_plink

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
fasindex="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta.fai"
REP="/PATH/04_fresh2/ref2.5repeats_header.bed"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/01_angsd_${1}
bgzip ${1}_${2}.vcf
tabix -p vcf ${1}_${2}.vcf.gz

# Extract biallelic regions from vcf containing variant and invariant sites
### generate variant sites only vcf ### NOT NEEDED ONLY VARIANTS CALLED
#bcftools view --threads 12 -V indels --min-alleles 2 --max-alleles 2 -Oz -o ${1}_$2_varonly.vcf.gz ${1}_${2}.vcf.gz
#bcftools index ${1}_${2}_varonly.vcf.gz

# Remove repeat regions
vcftools --exclude-bed ${REP} --gzvcf ${1}_${2}.vcf.gz --recode --out ${1}_${2}_norepeats
bcftools view -Oz -o ${1}_${2}_norepeats.vcf.gz ${1}_${2}_norepeats.recode.vcf
bcftools index ${1}_${2}_norepeats.vcf.gz
rm ${1}_${2}_norepeats.recode.vcf

#remove missingness
bcftools view --threads 12 -Oz -o ${1}_${2}_varonly_norepeats_miss20.vcf.gz --include \
         "F_MISSING <= 0" \
         ${1}_${2}_varonly_norepeats.vcf.gz
bcftools index ${1}_${2}_varonly_norepeats_nomiss.vcf.gz


#plink --vcf ${1}_${2}_varonly_norepeats_miss20.vcf.gz --allow-extra-chr --missing --out _norepeats_miss20


#then merge with outgroup
# run  2.2.1_outgroup.sh

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
