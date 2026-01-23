#!/bin/bash -l
#SBATCH -J blast
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 3.3_blast.sh D_RbY23
# sbatch 3.3_blast.sh DKoC53

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
module load blast-plus/2.9.0-gcc8

dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}

# cd ${dat}/03_assemblies/${1}_chr18_flye_ragtag
# makeblastdb -in ${dat}/ASM73873v6_chr18.fasta -dbtype nucl -out ASM73873v6_chr18
# blastn -query ${1}_chr18_flye_ragtag.fasta -db ASM73873v6_chr18 -out ${1}_blast_results.tsv \
# -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
# blastn -query ${1}_chr18_flye_ragtag.fasta -db ASM73873v6_chr18 -out ${1}_blast_pairwise2.txt -outfmt 0 -max_hsps 1

cd ${dat}/03_flye/${1}
#makeblastdb -in ${dat}/ASM73873v6_chr18.fasta -dbtype nucl -out ASM73873v6_chr18
#blastn -query ${1}_chr18_flye_inversion_region.fasta -db ASM73873v6_chr18 -out ${1}_inversion_region_blast_results.tsv \
#-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
#blastn -query ${1}_chr18_flye_invresion_region.fasta -db ASM73873v6_chr18 -out ${1}_inversion_region_blast_pairwise2.txt -outfmt 0 -max_hsps 1

blastn -query ${1}_chr18_flye_before_inversion_region.fasta -db ASM73873v6_chr18 -out ${1}_before_inversion_region_blast_results.tsv \
-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
blastn -query ${1}_chr18_flye_before_invresion_region.fasta -db ASM73873v6_chr18 -out ${1}_before_inversion_region_blast_pairwise2.txt -outfmt 0 -max_hsps 1

blastn -query ${1}_chr18_flye_after_inversion_region.fasta -db ASM73873v6_chr18 -out ${1}_after_inversion_region_blast_results.tsv \
-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
blastn -query ${1}_chr18_flye_after_invresion_region.fasta -db ASM73873v6_chr18 -out ${1}_after_inversion_region_blast_pairwise2.txt -outfmt 0 -max_hsps 1


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
