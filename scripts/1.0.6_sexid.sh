#!/bin/bash -l
#SBATCH -J sexid
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j.out

##sbatch 1.0.6_sexid.sh 00_eager_TE2
##sbatch 1.0.6_sexid.sh 00_eager_TE

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH/05_aDNA"

cd ${dat}/${1}
cd results/qualimap

for i in *_stats; do cd $i; grep "chr\|scaff" genome_results.txt | awk 'NR>2 {print $4}' >> ../01_cov_$i.txt; cd ../; done
for i in *_stats; do cd $i; grep "chr\|scaff" genome_results.txt | awk 'NR>2 {print $1}' > ../01_header.txt; done
for i in *stats; do base=${i%_rmdup_stats*}; echo ${base}; done > 01_samples.txt

paste 01_header.txt 01_cov_*.txt > 01_allchr_cov.txt
awk '                  
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' 01_allchr_cov.txt > 01_allchr_cov_transposed.txt

#sex-determination
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($386/$num<0.6) {print "female_" $386/$num}else {if ($386/$num>0.9) {print "male_" $386/$num} else {print "unsure_" $386/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff16_${i}.txt; done
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($442/$num<0.6) {print "female_" $442/$num}else {if ($442/$num>0.9) {print "male_" $442/$num} else {print "unsure_" $442/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff21_${i}.txt; done
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($563/$num<0.6) {print "female_" $563/$num}else {if ($563/$num>0.9) {print "male_" $563/$num} else {print "unsure_" $563/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff32_${i}.txt; done

paste 01_samples.txt sex_scaff*.txt > 01_sex.txt

