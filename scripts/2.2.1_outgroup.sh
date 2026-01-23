#!/bin/bash -l
#SBATCH -J outgroup
#SBATCH --cpus-per-task=12
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled outlier_geno_maxmis_q20_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_trans_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh wgs_rescaled all_geno_maxmis_q20_plink_varonly_norepeats moneduloides_monedula

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
out="/PATH/02_outgroups/07_call"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/01_angsd_${1}

bcftools query -f '%CHROM\t%POS\n' ${1}_${2}.vcf.gz > ${1}_${2}_positions.txt

#note the outgroup is arranged lexicographically but there is some positional sorting issue that bcftools did not fix
#outlier region is not affected as the wrongly order positions are not involved
#below id specific to neutral region where i did sort the extracted region of the outgroup vcf without AMcrow
# bcftools view --threads 12 -R ${1}_${2}_positions.txt ${out}/outgroups_noamcrow_incnonvariant.gvcf.gz | bgzip -c > outgroup_pos_selected_tmp.vcf.gz
# zcat outgroup_pos_selected_tmp.vcf.gz | grep '^#' > outgroup_pos_selected_tmp.vcf_header.txt
# zcat outgroup_pos_selected_tmp.vcf.gz | grep -v '^#' | sort -k1,1 -k2,2n > outgroup_pos_selected_tmp_sorted.vcf
# cat outgroup_pos_selected_tmp.vcf_header.txt outgroup_pos_selected_tmp_sorted.vcf | bgzip -c > outgroup_pos_selected_tmp_sorted.vcf.gz
# tabix -f -p vcf outgroup_pos_selected_tmp_sorted.vcf.gz
# bcftools merge --threads 12 ${1}_${2}.vcf.gz outgroup_pos_selected_tmp_sorted.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup.vcf.gz
# tabix -f -p vcf ${1}_${2}_${3}_outgroup.vcf.gz 
# rm outgroup_pos_selected_tmp*

#!!!this is used for admixtools!!!
bcftools view --threads 12 -V indels --min-alleles 2 --max-alleles 2 ${1}_${2}_${3}_outgroup.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
# bcftools view --threads 12 -Oz -o ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz --include "F_MISSING <= 20" ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
# tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz 

#for outlier/ wgs with no issue
bcftools view --threads 12 -R ${1}_${2}_positions.txt ${out}/outgroups_noamcrow_incnonvariant.gvcf.gz | bgzip -c > outgroup_pos_selected_tmp.vcf.gz
bcftools merge --threads 12 ${1}_${2}.vcf.gz outgroup_pos_selected_tmp.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup.vcf.gz
tabix -f -p vcf ${1}_${2}_${3}_outgroup.vcf.gz 
bcftools view --threads 12 -V indels --min-alleles 2 --max-alleles 2 ${1}_${2}_outgroup.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
bcftools view --threads 12 -Oz -o ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz --include "F_MISSING <= 20" ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz 

# for extraction of neutral sites with 3 outgroups, incamcrow (repeat from the first part)
# bcftools view --threads 12 -R ${1}_${2}_positions.txt ${out}/amcrow_moneduloides_monedula_incnonvariant.gvcf.gz | bgzip -c > outgroup_pos_selected_tmp.vcf.gz
# zcat outgroup_pos_selected_tmp.vcf.gz | grep '^#' > outgroup_pos_selected_tmp.vcf_header.txt
# zcat outgroup_pos_selected_tmp.vcf.gz | grep -v '^#' | sort -k1,1 -k2,2n > outgroup_pos_selected_tmp_sorted.vcf
# cat outgroup_pos_selected_tmp.vcf_header.txt outgroup_pos_selected_tmp_sorted.vcf | bgzip -c > outgroup_pos_selected_tmp_sorted.vcf.gz
# tabix -f -p vcf outgroup_pos_selected_tmp_sorted.vcf.gz
# bcftools merge --threads 12 ${1}_${2}.vcf.gz outgroup_pos_selected_tmp_sorted.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup.vcf.gz
# tabix -f -p vcf ${1}_${2}_${3}_outgroup.vcf.gz 
# rm outgroup_pos_selected_tmp*
# bcftools view --threads 12 -V indels --min-alleles 2 --max-alleles 2 ${1}_${2}_${3}_outgroup.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
# tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
# bcftools view --threads 12 -Oz -o ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz --include "F_MISSING <= 20" ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
# tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele_20miss.vcf.gz 


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
