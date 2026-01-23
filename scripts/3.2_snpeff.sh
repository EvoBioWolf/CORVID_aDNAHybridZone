#!/bin/bash -l
#SBATCH -J snpeff
#SBATCH --cpus-per-task=2
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 3.2_snpeff.sh

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}
cd $dat/02_results/uli_pca

#1: convert plink to vcf
plink --bfile add_fresh_anc_tmp_data_GG_crows --allow-extra-chr --recode vcf --chr-set 27 --chr 18 --from-bp 2090000 --to-bp 4591407 --out add_fresh_anc_tmp_data_GG_crows_peaks_7860
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {split($3, a, ":");$1 = a[1];$2 = a[2];print}' add_fresh_anc_tmp_data_GG_crows_peaks_7860.vcf | sed -e 's/^scaffold_78/NW_010959954.1/'| sed -e 's/^scaffold_60/NW_010959761.1/' > add_fresh_anc_tmp_data_GG_crows_peaks_modified_7860.vcf

conda activate java21_env
cd ${dat}/snpEff

#2: install data with NCBI genome: GCF_000738735.1
cd ./data
mkdir corvus_cornix
cd corvus_cornix
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_genomic.fna.gz
gunzip -c GCF_000738735.1_Hooded_Crow_genome_genomic.fna.gz > sequences.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_genomic.gff.gz
gunzip -c GCF_000738735.1_Hooded_Crow_genome_genomic.gff.gz > genes.gff #removed mitochondrial in the genes gff as it is problematic
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_protein.faa.gz
gunzip -c GCF_000738735.1_Hooded_Crow_genome_protein.faa.gz > protein.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_cds_from_genomic.fna.gz
gunzip -c GCF_000738735.1_Hooded_Crow_genome_cds_from_genomic.fna.gz > cds.fa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.1_Hooded_Crow_genome/GCF_000738735.1_Hooded_Crow_genome_rna.fna.gz
gunzip -c GCF_000738735.1_Hooded_Crow_genome_rna.fna.gz > rna.fa
rm*.gz

#3: build 
#modify config file to include corvus_cornix.genome
cd ${dat}/snpEff
java -Xmx4g -jar snpEff.jar build -gff3 -v corvus_cornix

#unfortunately I ran into issue with as the parentID of CDS (rna id) in the gff file is missing in the the header of cds.fa
#I will try to add rna id to the cds.fz file

grep CDS genes.gff | awk '{print $1,$3,$9}' | uniq > CDS_rna_ref_tmp.txt

awk '{
match($0, /Parent=([^;]+)/, parent); 
match($0, /Name=([^;]+)/, name); 
if (parent[1] != "" && name[1] != "") {
    print parent[1], name[1];
   }
}' CDS_rna_ref_tmp.txt | grep rna | uniq > CDS_rna_ref_table_tmp.txt

awk -v map=CDS_rna_ref_table_tmp.txt '
BEGIN {
  while ((getline < map) > 0) {
    rna_by_protein[$2] = $1;  # XP_* → rna*
  }
}
{
  if ($0 ~ /^>/) {
    match($0, /protein_id=([^]]+)/, p);
    rna = rna_by_protein[p[1]];
    if (rna != "") {
      print ">" rna "|" substr($0, 2);
    } else {
      print $0;
    }
  } else {
    print $0;
  } }' cds.fa > cds_with_rna.fa

awk '/^>/ {
  match($0, /^>(rna[0-9]+)/, m);
  if (m[1] != "") print ">" m[1];
  else print $0;  next}{ print }' cds_with_rna.fa > cds_clean.fa

mv cds.fa cds_unmodified.fs
mv cds_clean.fa cds.fa

cd ${dat}/snpEff
java -Xmx4g -jar snpEff.jar build -gff3 -v corvus_cornix
#if the database built succefully, a snpEffectPredictor.bin is produced

#4. annotate
java -Xmx4g -jar snpEff.jar -v corvus_cornix ${dat}/02_results/uli_pca/add_fresh_anc_tmp_data_GG_crows_peaks_modified_7860.vcf > ${dat}/02_results/uli_pca/add_fresh_anc_tmp_data_GG_crows_peaks_modified_7860_annotated.vcf

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/60)) mins to complete this task"
echo "It takes $((($ENDTIME - $STARTTIME)/3600)) hours to complete this task"

