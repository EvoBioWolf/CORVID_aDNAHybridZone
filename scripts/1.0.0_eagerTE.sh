#!/bin/bash -l
#SBATCH -J eager_TE
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.0_eagerTE.sh 00_eager_twist /PATH/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_twist
# sbatch 1.0.0_eagerTE.sh 00_eager_mybaits /PATH/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_mybaits

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30
module load nextflow
conda activate biotools

dat="/PATH/05_aDNA"

cd ${dat}
mkdir ${1}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input /PATH/05_aDNA/${3}.tsv --fasta ${2} \
-c /PATH/05_aDNA/base_modified2.config \
--snpcapture_bed /PATH/05_aDNA/probes_232015.bed \
--clip_adapters_list /PATH/05_aDNA/adapterlist.txt 

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
