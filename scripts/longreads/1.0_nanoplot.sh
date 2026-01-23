#!/bin/bash -l
#SBATCH -J nanoplot
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 1.0_nanoplot.sh 20250826_D_RbY23_UL D_RbY23_UL
# sbatch 1.0_nanoplot.sh 20250917_DNeY39_UL DNeY39_UL
# sbatch 1.0_nanoplot.sh 20250909_DKoC53_UL DKoC53_UL
# sbatch 1.0_nanoplot.sh 20251014_DLoC21_UL DKoC21_UL

echo $(date)
STARTTIME=$(date +%s)

conda activate nanoplot #changed, nanoplot removed
dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}
cd 00_rawreads

#tar -xvzf *.tar.gz
NanoPlot --fastq ${1}/fastq_pass/*.fastq.gz -o nanoplot_output_${2} --loglength --N50


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
