#!/bin/bash -l
#SBATCH -J recaleager
#SBATCH --cpus-per-task=8
#SBATCH --time=10-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.6_recal_eager.sh 00_eager_enrichedall /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_enrichedall_rescale.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_rescale.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_rescale_add.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall_02 /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_rescale_2.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall_02 /PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta /PATH/05_aDNA/acrow_eager_list_wgsall_rescale_add_2.tsv

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30

module purge
module load nextflow
conda activate biotools

dat="/PATH/05_aDNA"

cd ${dat}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input ${3} --fasta ${2} \
-c /PATH/05_aDNA/base_modified2.config \
--run_mapdamage_rescaling --skip_preseq --skip_deduplication

#merge the UDG and non-udg sequences post rescaling
# cd ./results/damage_rescaling
# samtools merge --threads 8 DVT014_UDG.trimmed_rescaled.bam DVT014_WGS.trimmed_rescaled.bam -o DVT014_merged.trimmed_rescaled.bam
# samtools merge --threads 8 WMP006_UDG.trimmed_rescaled.bam WMP006_WGS.trimmed_rescaled.bam -o WMP006_merged.trimmed_rescaled.bam
# samtools index -b --threads 8 DVT014_merged.trimmed_rescaled.bam
# samtools index -b --threads 8 WMP006_merged.trimmed_rescaled.bam


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
