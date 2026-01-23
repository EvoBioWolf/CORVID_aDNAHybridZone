#!/bin/bash -l
#SBATCH -J angsd1k
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.1_angsd_1111snp_AA.sh moneduloides 1111

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
out=" /PATH/02_outgroups/07_call"

cd ${dat}/02_results
cd ./uli_pca
mkdir 01_angsd_${1}
cd 01_angsd_${1}

vcftools --gzvcf ${out}/moneduloides_incnonvariant.gvcf.gz --positions ../knief_snp_1based_tabsep.pos --recode --out moneduloides_AA_1111
bgzip moneduloides_AA_1111.recode.vcf
tabix -p vcf moneduloides_AA_1111.recode.vcf.gz
bcftools query -f '%CHROM\t%POS\t[\t%TGT]\n' moneduloides_AA_1111.recode.vcf.gz | sed 's:/::g' |\
awk -v OFS="\t" 'BEGIN {print "CHROM", "POS", "FS66096", "NC3", "NC5", "NC6", "NC8"} {print}' >  moneduloides_AA_1111.txt

#hap_concensus_maxmis_q20
angsd -bam ${dat}/${3}.filelist -dohaplocall 2 -doCounts 1 -rf ../knief_snp_1based.pos \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 4 -minMinor 1 -out ${1}_${2}_hapconsensus_q20

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
