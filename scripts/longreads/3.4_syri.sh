#!/bin/bash -l
#SBATCH -J syri
#SBATCH --time=1-00:00:00
#SBATCH -o /PATH/06_longreads/slurms/slurm-%j-%x.out

# sbatch 3.4_syri.sh D_RbY23 D_RbY23_chr18_flye_ragtag
# sbatch 3.4_syri.sh DNeY39 DNeY39_chr18_flye_ragtag
# sbatch 3.4_syri.sh DKoC53 DKoC53_chr18_flye_ragtag
# sbatch 3.4_syri.sh DLoC21 DLoC21_chr18_flye_ragtag

echo $(date)
STARTTIME=$(date +%s)

conda activate syri_env

dat="/PATH/06_longreads"
ref="/PATH/06_longreads/GCF_000738735.6_ASM73873v6_genomic.mmi"

cd ${dat}

cd ${dat}/03_assemblies/${1}_chr18_flye_ragtag

#map to ref
# syri -c ${2}_ref.bam -r ${dat}/ASM73873v6_chr18.fna -q ${2}.fasta -k -F B --prefix ${1}_ref_
# plotsr --sr ${1}_ref_syri.out --genomes ./genomes.txt -o syri_${1}_ref.pdf -S 0.5 -W 10 -H 5 -f 8

#ref map to newly assembled chr18
# syri -c ref_to_${2}.bam -q ${dat}/ASM73873v6_chr18.fna -r ${2}.fasta -k -F B --prefix ref_to_${1}_
# plotsr --sr ref_to_${1}_syri.out --genomes ./genomes2.txt -o syri_ref_to_${1}.pdf -S 0.5 -W 10 -H 5 -f 8

#map to hooded RbY23 (3 alignments)
# syri -c ${1}_to_D_RbY23_chr18_flye_ragtag.bam -r ../D_RbY23_chr18_flye_ragtag/D_RbY23_chr18_flye_ragtag.fasta -q ${2}.fasta -k -F B --prefix ${1}_to_D_RbY23_
# plotsr --sr ../D_RbY23_chr18_flye_ragtag/D_RbY23_ref_syri.out --sr ${1}_to_D_RbY23_syri.out --genomes ./genomes3.txt -o syri_ref_D_RbY23_${1}.pdf -S 0.5 -W 10 -H 5 -f 8

#map to RbY23 (2 alignments)
# plotsr --sr ${1}_to_D_RbY23_syri.out --genomes ./genomes2_2.txt -o syri_D_RbY23_${1}.pdf -S 0.5 -W 10 -H 5 -f 8

#map to DNeY39 (2 alignments)
# syri -c ${1}_to_DNeY39_chr18_flye_ragtag.bam -r ../DNeY39_chr18_flye_ragtag/DNeY39_chr18_flye_ragtag.fasta -q ${2}.fasta -k -F B --prefix ${1}_to_DNeY39_
# plotsr --sr ${1}_to_DNeY39_syri.out --genomes ./genomes2_3.txt -o syri_DNeY39_${1}.pdf -S 0.5 -W 10 -H 5 -f 8

# map to DKoC53
# syri -c ${1}_to_DKoC53_chr18_flye_ragtag.bam -r ../DKoC53_chr18_flye_ragtag/DKoC53_chr18_flye_ragtag.fasta -q ${2}.fasta -k -F B --prefix ${1}_to_DKoC53_
# plotsr --sr ${1}_to_DKoC53_syri.out --genomes ./genomes2_4.txt -o syri_DKoC53_${1}.pdf -S 0.5 -W 10 -H 5 -f 8

#alignment series up to 4 in order: oriRef > D_RbY23 > DNeY39 > DKoC53 > DLoc21
#final plot used
syri -c ../DKoC53_chr18_flye_ragtag/DKoC53_to_DLoc21_chr18_flye_ragtag.bam -r ${2}.fasta -q ../DKoC53_chr18_flye_ragtag/DKoC53_chr18_flye_ragtag.fasta -k -F B --prefix DKoC53_to_DLoc21_
plotsr --sr ../D_RbY23_chr18_flye_ragtag/D_RbY23_ref_syri.out --sr ../DNeY39_chr18_flye_ragtag/DNeY39_to_D_RbY23_syri.out \
--sr DLoc21_to_DNeY39_syri.out --sr DKoC53_to_DLoc21_syri.out \
--genomes ./genomes4_2.txt -o syri_D_RbY23_DNeY39_DLoc21_DKoC53.pdf -S 0.5 -W 10 -H 5 -f 8


ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
