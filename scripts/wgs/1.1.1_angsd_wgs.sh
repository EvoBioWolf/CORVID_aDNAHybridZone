#!/bin/bash -l
#SBATCH -J WGSangsd
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

## sbatch 1.1.1_angsd_wgs.sh wgs_rescaled all
## sbatch 1.1.1_angsd_wgs.sh wgs_rescaled_outgroup all

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
mkdir 01_angsd_${1}
cd ${dat}/01_angsd_${1}

missIND=7
minIND=71
#wgsallfresh_rescaled: 130 ind: 17 ancient + 113 fresh, dp3, missIND=25 - over 90M SNPs
#wgsrescaledoutrgoup: 78 ind: 13 anc + 65 fresh, dp=3, missIND=7 if missing in more than 50% of the anc samples then remove sites, maf=2

#hap_concensus_maxmis_q20
#angsd -bam bam_${1}.filelist -dohaplocall 2 -doCounts 1 \
#-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
#-nThreads 24 -minMinor 2 -maxMis ${missIND} -setMinDepthInd 3 -out ${1}_${2}_hapconsensus_maxmis_q20_dp3

#hap_concensus_maxmis_q20
# angsd -bam bam_${1}.filelist -dohaplocall 2 -doCounts 1 \
# -checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
# -nThreads 16 -minMinor 2 -maxMis ${missIND} -out ${1}_${2}_hapconsensus_maxmis_q20

#geno
#angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
#-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -minInd ${minIND} -remove_bads 1 \
#-nThreads 24 -setMinDepthInd 3 -geno_minDepth 3 -out ${1}_${2}_geno_maxmis_q20_dp3

#geno_chr18
-rf ../knief_snp_1based.pos


#zless wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20.vcf.gz |\
#awk '$1 != "scaffold_60" && $1 != "scaffold_78"' > wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_NOoutlier.vcf
#gzip wgs_rescaled_outgroup_all_hapconsensus_maxmis_q20_NOoutlier.vcf

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
