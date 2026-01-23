#!/bin/bash -l
#SBATCH -J mapref5
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# cd /PATH/04_fresh2/04_markdup
# for i in C.cornix_IRQ*_markdup_cleaned.ploidy2.bam; do base=${i%_markdup_cleaned.ploidy2.bam*}; sbatch /PATH/05_aDNA/4.0_ref5.7_fresh.sh ${base}; done #notused

# for i in C.corone_E*_markdup_cleaned.ploidy2.bam; do base=${i%_markdup_cleaned.ploidy2.bam*}; sbatch /PATH/05_aDNA/4.0_ref5.7_fresh.sh ${base}; done
# for i in C.cornix_B*_markdup_cleaned.ploidy2.bam; do base=${i%_markdup_cleaned.ploidy2.bam*}; sbatch /PATH/05_aDNA/4.0_ref5.7_fresh.sh ${base}; done
# for i in C.corone_D*_markdup_cleaned.ploidy2.bam; do base=${i%_markdup_cleaned.ploidy2.bam*}; sbatch /PATH/05_aDNA/4.0_ref5.7_fresh.sh ${base}; done

# for i in C.cornix_IRQ*_markdup_cleaned.ploidy2.bam; do base=${i%_markdup_cleaned.ploidy2.bam*}; sbatch /PATH/05_aDNA/4.0_ref5.7_fresh.sh ${base}; done

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/05_aDNA"
ref="/PATH/06_longreads/ASM73873v6"

cd $dat
cd 00_eager_enrichedall/results
cd mapped_ref5.7

for i in {1..15} {17..24} {26..28} 1A 4A
do
echo "Processing chr${i}"
samtools view -b -@ 8 /PATH/04_fresh2/04_markdup/${1}_markdup_cleaned.ploidy2.bam \
$(cat /PATH/01_probes/scaffold2chr/chr${i}.scaffolds) > ${1}_chr${i}_markdup_cleaned.ploidy2.bam
samtools index ${1}_chr${i}_markdup_cleaned.ploidy2.bam
samtools fastq -F 0x400 -@ 8 ${1}_chr${i}_markdup_cleaned.ploidy2.bam | pigz -p 8 > ${1}_chr${i}.fastq.gz #exclude duplicate reads
bwa mem -t 8 -M ${ref}_chr${i}.fasta ${1}_chr${i}.fastq.gz | samtools view -bhq 6 -@ 8 | samtools sort -m 4G -@ 8 -o ${1}_ref5.7_chr${i}.bam
samtools index ${1}_ref5.7_chr${i}.bam
done

samtools merge -@ 8 ${1}_all_chr_ref5.7.bam ${1}_ref5.7_chr*.bam
samtools index ${1}_all_chr_ref5.7.bam

#rename after all done
#rename C.corone_E E C.corone_E*
#rename C.cornix_IRQ IRQ C.cornix_IRQ*
#rename C.cornix_B B C.cornix_B*
#rename C.corone_D D C.corone_D*

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

