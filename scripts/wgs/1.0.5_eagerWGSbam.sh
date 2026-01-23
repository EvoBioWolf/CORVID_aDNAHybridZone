#!/bin/bash -l
#SBATCH -J eagerWGSbam
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.5_eagerWGSbam.sh 00_eager_wgsall /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_bam.tsv
# sbatch 1.0.5_eagerWGSbam.sh 00_eager_wgsall_02 /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_bam.tsv

echo $(date)
STARTTIME=$(date +%s)

#module load openjdk 
#module load flux-core
module load charliecloud/0.30
#module purge
module load nextflow
conda activate biotools

dat="/PATH/05_aDNA"
projectDir="/PATH/05_aDNA/00_eager_wgsall_02"
cd ${dat}

#need to make sure bam file names are unique
# cd ${dat}/00_eager_wgs3/results/ext_mapped
# for i in *_sorted_noqc.bam
# do
# base=${i%_sorted_noqc.bam*}
# mv $i ${base}_3_sorted_noqc.bam 
# done

# cd ${dat}/00_eager_wgs2/results/ext_mapped 
# for i in *_sorted_noqc.bam
# do
# base=${i%_sorted_noqc.bam*}
# mv $i ${base}_2_sorted_noqc.bam
# done

# cd ${dat}/00_eager_wgs1/results/ext_mapped
# for i in *_sorted_noqc.bam
# do
# base=${i%_sorted_noqc.bam*}
# mv $i ${base}_1_sorted_noqc.bam
# done

cd ${dat}
mkdir ${1}
cd $1

unset DISPLAY
nextflow run nf-core/eager -r 2.4.7 -profile conda --input ${3} --fasta ${2} \
-c /PATH/05_aDNA/base_modified2.config \
--skip_damage_calculation \
--run_trim_bam --bamutils_clip_single_stranded_half_udg_left 1 --bamutils_clip_single_stranded_half_udg_right 3 \
--bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5 \
--run_genotyping --genotyping_tool angsd \
--write_allele_frequencies --angsd_glformat beagle --mtnucratio_header chrM --run_mtnucratio

#rerun this again as i realised half udg treated samples were not trimmed

#turned off mapdamage_rescaling and damage_calculation
ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
