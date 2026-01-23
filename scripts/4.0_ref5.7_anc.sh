B#!/bin/bash -l
#SBATCH -J mapref5
#SBATCH --get-user-env
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out


#sbatch 4.0_ref5.7_anc.sh HOC005
#sbatch 4.0_ref5.7_anc.sh ENG001
#sbatch 4.0_ref5.7_anc.sh WMP005
#sbatch 4.0_ref5.7_anc.sh WMP006
#sbatch 4.0_ref5.7_anc.sh RTT001
#sbatch 4.0_ref5.7_anc.sh GAB001 
#sbatch 4.0_ref5.7_anc.sh LSS004

#sbatch 4.0_ref5.7_anc.sh DVT014
#sbatch 4.0_ref5.7_anc.sh DVT016
#sbatch 4.0_ref5.7_anc.sh DVT022
#sbatch 4.0_ref5.7_anc.sh KZR002
#sbatch 4.0_ref5.7_anc.sh NCP001
#sbatch 4.0_ref5.7_anc.sh NCP002

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/05_aDNA"
ref="/PATH/06_longreads/ASM73873v6"

# I need to produce one ref per chr
# sed -E 's/>NC_[0-9]+.*chromosome ([0-9]+[A-Z]*).*/>chr\1/' GCF_000738735.6_ASM73873v6_genomic.fna > GCF_000738735.6_ASM73873v6_genomic_renamed.fna
# awk '/^>/{ if (seq) {print seq > out; } 
# seq=""; chr=$1; sub(/^>/,"",chr); out="ASM73873v6_"chr".fasta"; 
# print $0 > out; next} { print $0 >> out }' GCF_000738735.6_ASM73873v6_genomic_renamed.fna
# for i in ASM73873v6_chr*.fasta; do samtools faidx $i; done
# for i in ASM73873v6_chr*.fasta; do bwa index $i; done 

cd $dat
cd 00_eager_enrichedall/results
mkdir mapped_ref5.7
cd mapped_ref5.7

for i in {1..15} {17..24} {26..28} 1A 4A
do
echo "Processing chr${i}"
samtools view -b -@ 8 ../damage_rescaling/${1}_libmerged.trimmed_rescaled.bam $(cat /PATH/01_probes/scaffold2chr/chr${i}.scaffolds) > ${1}_chr${i}_libmerged.trimmed_rescaled.bam
samtools index ${1}_chr${i}_libmerged.trimmed_rescaled.bam
samtools fastq -F 0x400 -@ 8 ${1}_chr${i}_libmerged.trimmed_rescaled.bam | pigz -p 8 > ${1}_chr${i}.fastq.gz #exclude duplicate reads
bwa mem -t 8 -k 15 -B 2 -M ${ref}_chr${i}.fasta ${1}_chr${i}.fastq.gz | samtools view -bhq 6 -@ 8 | samtools sort -m 4G -@ 8 -o ${1}_ref5.7_chr${i}.bam
samtools index ${1}_ref5.7_chr${i}.bam
done

samtools merge -@ 8 ${1}_all_chr_ref5.7.bam ${1}_ref5.7_chr*.bam
samtools index ${1}_all_chr_ref5.7.bam

#samtools merge -@ 8 ${1}_all_chr_ref5.7.bam \
#${1}_chr1.bam ${1}_chr2.bam ${1}_chr3.bam ${1}_chr4.bam ${1}_chr5.bam ${1}_chr6.bam ${1}_chr7.bam ${1}_chr8.bam ${1}_chr9.bam \
#${1}_chr10.bam ${1}_chr11.bam ${1}_chr12.bam ${1}_chr13.bam ${1}_chr14.bam ${1}_chr15.bam ${1}_chr17.bam ${1}_chr18.bam ${1}_chr19.bam ${1}_chr20.bam \
#${1}_chr21.bam ${1}_chr22.bam ${1}_chr23.bam ${1}_chr24.bam ${1}_chr26.bam ${1}_chr27.bam ${1}_chr28.bam ${1}_chr1A.bam ${1}_chr4A.bam

#bam aln does not work for remapping to the new ref due to too short reads 

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"

