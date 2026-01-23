#!/bin/bash -l
#SBATCH -J subsetsnp
#SBATCH --cpus-per-task=6
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.3_subsetvariants.sh enrichedallfresh_rescaled probes_segAM72k_neutral_1based enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink_filtered_norm_moneduloides_monedula_outgroup_biallele segAM72k
# sbatch 1.3_subsetvariants.sh enrichedallfresh_rescaled probes_segAM72k_neutral_1based enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink_filtered segAM72k
# sbatch 1.3_subsetvariants.sh enrichedallfresh_rescaled_outgroup probes_segAM72k_neutral_1based neutral_hapconsensus_maxmis_q20 neutral_seg72k_hapconsensus_maxmis_q20
# sbatch 1.3_subsetvariants.sh enrichedallfresh_rescaled_outgroup probes_segAM72k_neutral_1based neutral_hapconsensus_maxmis_q20_dp3 neutral_seg72k_hapconsensus_maxmis_q20_dp3

# sbatch 1.3_subsetvariants.sh enrichedallfresh_rescaled_outgroupAM probes_segAM72k_neutral_1based neutral_hapconsensus_maxmis_q20 neutral_seg72k_hapconsensus_maxmis_q20

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
fasindex="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta.fai"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/01_angsd_${1}
vcftools --positions ${dat}/${2}.pos --gzvcf ${1}_${3}.vcf.gz --recode --out ${1}_${4}
mv ${1}_${4}.recode.vcf ${1}_${4}.vcf
bgzip ${1}_${4}.vcf
tabix -p vcf ${1}_${4}.vcf.gz

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"
