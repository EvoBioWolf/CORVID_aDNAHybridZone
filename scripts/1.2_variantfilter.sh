#!/bin/bash -l
#SBATCH -J filter
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.2_variantfilter.sh enrichedallfresh_rescaled_filtered_rus neutral_geno_maxmis_q20_dp3_plink
# sbatch 1.2_variantfilter.sh enrichedallfresh_rescaled_filtered_rus outlier_geno_maxmis_q20_dp3_plink
# sbatch 1.2_variantfilter.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_plink
# sbatch 1.2_variantfilter.sh enrichedallfresh_rescaled outlier_geno_maxmis_q20_plink
# sbatch 1.2_variantfilter.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_trans_plink

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
fasindex="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta.fai"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/01_angsd_${1}

# bgzip ${1}_${2}.vcf
# tabix -p vcf ${1}_${2}.vcf.gz

#remove the problematic samples
bcftools view --samples ^BDG002_TE_GB,DSZ007_TE_PL,DVT017_TE_B,KCZ001_TE_PL ${1}_${2}.vcf.gz | bgzip -c > ${1}_${2}_filtered.vcf.gz 
tabix -p vcf ${1}_${2}_filtered.vcf.gz

# Sort the vcf file for neutral angsd output (outlier is already sorted lexicographically)
mv ${1}_${2}_filtered.vcf.gz ${1}_${2}_filtered_unsorted.vcf.gz
mv ${1}_${2}_filtered.vcf.gz.tbi ${1}_${2}_filtered_unsorted.vcf.gz.tbi
bcftools view ${1}_${2}_filtered_unsorted.vcf.gz | grep -v "^#" | sort -k1,1 -k2,2n | bgzip -c > ${1}_${2}_filtered_sorted_tmp.vcf.gz
zcat ${1}_${2}_filtered_sorted_tmp.vcf.gz | cut -f1 | uniq > ${1}_${2}_sorted_scaffolds_tmp.txt
awk 'NR==FNR {len[$1]=$2; next} {if ($1 in len) print "##contig=<ID=" $1 ",length=" len[$1] ">"}' ${fasindex} ${1}_${2}_sorted_scaffolds_tmp.txt > ${1}_${2}_new_contigs_tmp.txt
zgrep -v "^##contig=" ${1}_${2}_filtered_unsorted.vcf.gz | grep "^#" > ${1}_${2}_header_without_contigs_tmp.txt
{
    head -n 4 ${1}_${2}_header_without_contigs_tmp.txt
    cat ${1}_${2}_new_contigs_tmp.txt
    tail -n +5 ${1}_${2}_header_without_contigs_tmp.txt 
} > ${1}_${2}_updated_header_tmp.txt
cat ${1}_${2}_updated_header_tmp.txt <(zcat ${1}_${2}_filtered_sorted_tmp.vcf.gz) | bgzip -c > ${1}_${2}_filtered.vcf.gz
tabix -f -p vcf ${1}_${2}_filtered.vcf.gz
rm ${1}_${2}_filtered_sorted_tmp.vcf.gz

#correct ref allele (angsd is not reliable)
bcftools norm --check-ref s -f ${ref} ${1}_${2}_filtered.vcf.gz | bgzip -c > ${1}_${2}_filtered_norm.vcf.gz
tabix -f -p vcf ${1}_${2}_filtered_norm.vcf.gz
#REF/ALT total/modified/added:   36477/4094/6 (outlier)
#REF/ALT total/modified/added:   191085/16233/19 (neutral)


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
