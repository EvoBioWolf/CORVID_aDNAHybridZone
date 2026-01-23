#!/bin/bash -l
#SBATCH -J eagerTEbam
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out


# sbatch 1.0.5_eagerTEbam.sh 00_eager_enrichedall /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_enrichedall_bam.tsv probes_232015_SNPsite

# sbatch 1.0.5_eagerTEbam.sh 00_eager_mybaits /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_mybaits_bam.tsv probes_104k

echo $(date)
STARTTIME=$(date +%s)

module load purge
module load charliecloud
module load nextflow/21.04.0-gcc8
conda activate biotools

dat="/PATH/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

#with snpsite bed files and trimmed bam
unset DISPLAY
# nextflow run nf-core/eager -profile conda -r 2.4.7 --input ${3} --fasta ${2} \
# -c /PATH/05_aDNA/base_modified2.config \
# --snpcapture_bed /PATH/05_aDNA/${4}.bed \
# --mtnucratio_header chrM --run_mtnucratio \
# --run_trim_bam --bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5 

# include duplicates
nextflow run nf-core/eager -profile conda -r 2.4.7 --input ${3} --fasta ${2} \
-c /PATH/05_aDNA/base_modified2.config \
--snpcapture_bed /PATH/05_aDNA/${4}.bed \
--skip_preseq --skip_deduplication --skip_damage_calculation \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5

#changed to a smaller probe set specific to mybaits (probes_104k.bed) to account coverage and depth of mybaits seq
# turned off mapdamage_rescaling for enriched_all and bam mapping quality threshold was introduced in enriched_all 
# --run_mapdamage_rescaling

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
