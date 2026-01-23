#!/bin/bash -l
#SBATCH -J map
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 1.1_map.sh 2025AUG26_D_RbY23_UL/fastq_pass D_RbY23 
# sbatch 1.1_map.sh 20250917_DNeY39_UL/fastq_pass DNeY39
# sbatch 1.1_map.sh 20250909_DKoC53_UL/fastq_pass DKoC53
# sbatch 1.1_map.sh 20251014_DLoC21_UL/fastq_pass DKoC21

echo $(date)
STARTTIME=$(date +%s)


dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}
cd 01_mapped

minimap2 -d ${dat}/GCF_000738735.6_ASM73873v6_genomic.mmi ${dat}/GCF_000738735.6_ASM73873v6_genomic.fna

minimap2 -ax map-ont -t 8 ${ref} ${dat}/00_rawreads/${1}/*.fastq.gz | samtools view -b -@ 8 | samtools sort -@ 8 -o ${2}.sorted.bam

samtools view -b -f 4 -@ 8 ${2}.sorted.bam > ${2}.unmapped.bam

qualimap bamqc -bam ${2}.sorted.bam -c --java-mem-size=32G

conda activate nanoplot
NanoPlot --bam ${2}.sorted.bam -o nanoplot_${2}_bam_output --N50 --loglength --title "Mapped Reads QC" --plots hex dot kde

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
