#!/bin/bash -l
#SBATCH -J COMangsd
#SBATCH --cpus-per-task=12
#SBATCH --time=10-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled neutral probes_193768_neutral.1based
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled neutral_seg72k probes_segAM72k_neutral.1based
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled outlier_1111 knief_snp_1based.pos
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled outlier probes_38498_outlier.1based

#outgroup for admixtool analysis
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled_outgroup neutral probes_193768_neutral.1based #haponly with dumpcount
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled_outgroup outlier_1111 knief_snp_1based.pos #haponly dumpcount


echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
# awk '{print $1":"$2}' /PATH/01_probes/snp_panel_2/final_snp_panel/neuall_outgroup_ascertained_nooverlap_final_TV_hwe.txt > probes_segAM72k_neutral.1based

cd ${dat}
mkdir 01_angsd_${1}

if [ $(ls ./01_angsd_${1}/bam_*.filelist | wc -l) -eq 1 ]
then
    echo "bam filelist already exist"
else
    echo "generate bam filelist"
    cp bam_${1}.filelist ./01_angsd_${1}/bam_${1}.filelist
fi

cd ${dat}/01_angsd_${1}

missIND=25
minIND=139
#enrichedallfresh_rescaled: 164 ind: 51 ancient + 113 fresh, no dp filtered, and dp3
#dp3 cause a lot of missingness in anc but no dp cause excussive alt sites

angsd -bam bam_${1}.filelist -dohaplocall 2 -doCounts 1 -rf ${dat}/${3} \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 12 -minMinor 1 -maxMis ${missIND} -out ${1}_${2}_hapconsensus_maxmis_q20 -dumpCounts 2
# then run haploidconversion.ipynb to convert to vcf 

#geno
angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -minInd ${minIND} -remove_bads 1 -geno_minDepth 3 \
-nThreads 12 -rf ${dat}/${3} -out ${1}_${2}_geno_maxmis_q20_dp3

mv ${1}_${2}_geno_maxmis_q20_dp3.tfam ${1}_${2}_geno_maxmis_q20_dp3_original.tfam
awk '{print $2, $1, "0", "0", "0", "-9"}' pop.txt >> ${1}_${2}_geno_maxmis_q20_dp3.tfam
plink --tfile ${1}_${2}_geno_maxmis_q20_dp3 --allow-no-sex --allow-extra-chr -make-bed --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --recode vcf --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --pca --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --missing --out ${1}_${2}_geno_maxmis_q20_dp3_plink
bgzip ${1}_${2}_geno_maxmis_q20_dp3_plink.vcf
tabix -p vcf ${1}_${2}_geno_maxmis_q20_dp3_plink.vcf.gz


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"Y
