#!/bin/bash -l
#SBATCH -J admixtools
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_plink_filtered_norm_moneduloides_monedula_outgroup_biallele enriched_neutral_160inds 1 41 50 8pop
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_trans_plink_filtered_norm_moneduloides_monedula_outgroup_biallele enriched_neutral_trans_160inds 1 1 20 9pop
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_plink_filtered_norm_moneduloides_monedula_outgroup_biallele_seg72k.recode enriched_segAM72k_160inds 1 41 50 8pop

# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup outlier_1111_hapconsensus_maxmis_q20_dp3 outlier_dp3
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup outlier_1111_hapconsensus_maxmis_q20 outlier

# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup neutral_hapconsensus_maxmis_q20_dp3 neutral_all_dp3
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup neutral_seg72k_hapconsensus_maxmis_q20_dp3 neutral_seg72k_dp3

# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup neutral_hapconsensus_maxmis_q20 neutral_all             
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroup neutral_seg72k_hapconsensus_maxmis_q20 neutral_seg72k

# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroupAM neutral_hapconsensus_maxmis_q20 neutral_all_AM
# sbatch 2.2.2_admixtools.sh enrichedallfresh_rescaled_outgroupAM neutral_seg72k_hapconsensus_maxmis_q20 neutral_seg72k_AM

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

conda activate biotools
echo $(date)
STARTTIME=$(date +%s)

# I first need to sort the vcf file
cd ${dat}/01_angsd_${1}
zless ${1}_${2}.vcf.gz | grep "^#" > ${1}_${2}.header.txt
zless ${1}_${2}.vcf.gz | grep -v "^#" | sort -k1,1 -k2,2n > ${1}_${2}.body.txt
cat ${1}_${2}.header.txt ${1}_${2}.body.txt | bgzip -c > ${1}_${2}.sorted.vcf.gz
rm ${1}_${2}.header.txt ${1}_${2}.body.txt

conda activate base #python compatibility
cd ${dat}/02_results/admixtools
python vcf2eigenstrat.py -v ${dat}/01_angsd_${1}/${1}_${2}.sorted.vcf.gz -o $3

cp ${3}.geno ${3}_pop.geno

#cp enriched_neutral_outgroupAM_pop.ind ${3}_pop.ind
cp enriched_neutral_outgroup_pop.ind ${3}_pop.ind

awk '{print $2":"$4, $2, $3, $4, $5, $6}' ${3}.snp | awk '{gsub("scaffold_", "", $2)}1' | awk '{print $1, $2+1, $4*0.000001, $4, $5, $6}' OFS="\t" > ${3}_final.snp
mv ${3}_final.snp ${3}_pop.snp

# conda activate r4.2
# cd ${dat}/02_results/admixtools
# Rscript ./admixtools_1000gen.R ${3}_pop $4 $5 $6 $7

#cat *bestscores.txt | grep -v NA | grep -v run_no | sort -n -k2,2g

#summary & compile
# mkdir 01_${2}
# mv ${2}*_adm* ./01_${2}
# cd 01_${2} 
# for i in *_bestscores.txt; do echo $i >> 01_summary.txt; awk 'NR>1 {print $2, $1}' $i | sort -n | head -1 >> 01_summary.txt; done
# mkdir 01_bestgraphs
# for i in *_bestscores.txt; do base=${i%_bs1to100_bestscores.txt*}; sort -k2 -n $i | awk 'NR>1 {print $1}' | cp $(echo ${base}_run_$(head -n 1).pdf) ./01_bestgraphs; done
# conda activate base
# cd 01_bestgraphs
# for i in *.pdf; do base=${i%_adm*}; convert -append ${base}_*.pdf merged_${base}.pdf; done

