# Pipeline for ancient crow project
All scripts are currently stored in the github [folder](./scripts/). Rscripts for plotting are stored [here](./scripts/plot/).

# nfcore/Eager2

I ran [Eager](https://nf-co.re/eager/2.4.7) version 2.4.7 on BioHPC wih conda profile. The Eager pipeline works with nextflow and is designed for ancient DNA analysis.
I used the pipeline for target-enriched data (as well as WGS for checking). For target-enriched data a SNP bed file was provided in addition to the reference genome v2.5 with chrW.

I broke the Eager pipeline into two parts. First raw reads are processed up to the part just before mapping. I mapped the trimmed reads to the reference genome with my own script using BWA-aln to multithread more efficiently. I then resume the Eager pipeline with bam input (sorted but duplicates have not been removed).

---

## Part 1: raw reads processing

Run [1.0.0_eagerTE.sh](./scripts/1.0.0_eagerTE.sh) for target-enriched
```
#!/bin/bash -l
#SBATCH -J eager_TE
#SBATCH --cpus-per-task=8
#SBATCH --time=14-00:00:00
#SBATCH -o PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 1.0.0_eagerTE.sh 00_eager_twist /PATH/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_twist
# sbatch 1.0.0_eagerTE.sh 00_eager_mybaits /PATH/04_fresh2/genome_HC_allpaths41687_v2.5_chrW.fasta acrow_eager_list_mybaits

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30
module load nextflow
conda activate biotools

dat="/PATH"

cd ${dat}
mkdir ${1}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input /PATH/${3}.tsv --fasta ${2} \
-c /PATH/base_modified2.config \
--snpcapture_bed /PATH/probes_232015.bed \
--clip_adapters_list /PATH/adapterlist.txt

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
```

An input file *.tsv is required to run the pipeline on multiple samples. It has a header like below. Check the [Eager page](https://nf-co.re/eager/2.5.1/docs/usage) for more info:

```
 Sample_Name    Library_ID  Lane    Colour_Chemistry    SeqType Organism    Strandedness    UDG_Treatment   R1  R2  BAM
```

[adapterlist.txt](./scripts/adapterlist.txt) is provided to specify 7 pairs of adapter combinations to be trimmed off from R1 and R2. The read 1 adapter used to make our [ssDNA libraries](https://www.nature.com/articles/s41596-020-0338-0) is 5bp shorter than the default Illumina adapter, so we need to provide this unique adapter pair to AdapterRemoval (shown in the 1st and 2nd rows, while the 3rd and 4th rows show the default pair of adapters). We also included the 6th and 7th pairs of seqeunces as Fastqc flagged over-representation of these unkown sequences in multiple samples after trimming. I have blasted these sequences and they seem to be bacteria or cloning vector.

```
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  GGAAGAGCGTCGTGTAGGGAAAGAGTGT
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG  GGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNAGATCTCGGTGGTCGCCGTATCATT
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNAGATCTCGGTGGTCGCCGTATCATT
AGATCGGAAGA GGAAGAGCGTCG
ATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  TCTTCCGATCTGGAAGAGCGTCGTGTAGGGAAAGAGTGT
ATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGA    TCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTG
```

AdapterRemoval will also merge read1 and read2 together to form longer reads. The merged reads are found for each sequenced sample in the folder `./results/adapterremoval/output`. Samples of the same libraries, but sequenced on different sequencing lane(s) will be merged and are found in the folder `./results/lanemerging`

I stopped the pipeline after AdapterRemoval and lanemerging (for the wgs samples).

---

## Part 2: BWA-aln mapping

Run [1.0.1_mapTE.sh](./scripts/1.0.1_mapTE.sh) for target-enriched data. Seeding is disabled to [optimized mapped reads](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3468387/), but more computation time is required. This is the main reason why I stopped the Eager pipeline to map each sample as a separate job - increasing the number of CPUs on nextflow config file did not improve mapping efficiency

[1.0.1_mapTE.sh](./scripts/1.0.1_mapTE.sh)

```
#!/bin/bash -l
#SBATCH -J mapTE
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/slurms/slurm-%j-%x.out

# cd adapterremoval/output
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $1"_"$2}'); do sbatch /PATH/1.0.1_mapTE.sh $i 00_eager_twist adapterremoval/output; done
# for i in $(ls *.fq.gz | awk -v FS="_" '{print $2"_"$3}'); do sbatch /PATH/1.0.1_mapTE.sh $i 00_eager_mybaits adapterremoval/output; done

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="/PATH"
ref="/PATH/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd ${2}/results
mkdir ext_mapped
cd ${3}

echo ${1}

for fname in *${1}*.fq.gz
do
bwa aln -l 1024 -n 0.04 -o 2 -t 8 ${ref} ${fname} > ${dat}/${2}/results/ext_mapped/${1}.sai
done
#default seeding for -l is 32bp   
#-n changed from 0.02 to 0.04 (stricter)

cd ${dat}/${2}/results/ext_mapped
for fname in ${1}.sai
do
bwa samse ${ref} ${1}.sai ${dat}/${2}/results/${3}/*${1}*.fq.gz > ${1}.sam
done

#no quality threshold
samtools view -bh -@ 8 -bS ${1}.sam > ${1}_noqc.bam
samtools sort -@ 8 ${1}_noqc.bam -o ${1}_sorted_noqc.bam
rm ${1}_noqc.bam
```
`1.0.1_mapWGS.sh`, which is the same as above, used for WGS data

---

## Part 3: BAM processing

Run [1.0.5_eagerTEbam.sh](./scripts/1.0.5_eagerTEbam.sh) for target-enriched data.

Continuing the Eager pipeline, but with sorted BAM files instead of fastq files this time.

[1.0.5_eagerTEbam.sh](./scripts/1.0.5_eagerTEbam.sh)

```
#!/bin/bash -l
#SBATCH -J eagerTEbam
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o /PATH/slurms/slurm-%j-%x.out

# sbatch 1.0.5_eagerTEbam.sh 00_eager_enrichedall /PATH/genome_HC_allpaths41687_v2.5.fasta /PATH/acrow_eager_list_enrichedall_bam.tsv

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30
module purge
module load nextflow
conda activate biotools

dat="/PATH"

cd ${dat}
mkdir ${1}
cd $1

#with snpsite bed files and trimmed bam
unset DISPLAY

nextflow run nf-core/eager -r 2.4.7 -profile conda --input ${3} --fasta ${2} \
-c /PATH/base_modified2.config \
--snpcapture_bed /PATH/probes_232015_SNPsite.bed \
--mtnucratio_header chrM --run_mtnucratio \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 5 --bamutils_clip_single_stranded_none_udg_right 5 \
--run_genotyping --genotyping_tool angsd --pileupcaller_bedfile /PATH/probes_232015_SNPsite.bed \
--write_allele_frequencies --angsd_glformat beagle

# turned off mapdamage_rescaling for enriched_all and bam mapping quality threshold was introduced in enriched_all
# --run_mapdamage_rescaling
```


Like the first part of Eager pipeline, we need a *.tsv input file to specify where the bam files are located. All bam files need to have a unique name even if the paths are different.

The pipeline runs multiple steps and it is important to know the sequence of events after the sorted bam files have been provided. Below is the chronological order of folders being produced at each step:

1. Index input bam
2. Samtools flagstat: % of mapped reads
3. Preseq: duplication rate
4. Endorspy: percentage of endogenous DNA
5. Markduplicates: removes duplicates
6. Mtnucratio: mtDNA to nuclear ratio
7. Damageprofiler: damage calculation (this step often runs into problem for wgs samples)
8. Seqtype_merge: merge bam files of same sample and library type (partial or no UDG seperately) (if needed)
9.  Mapdamage_rescaling: probablistically replaces Ts back to Cs based on reference-mismatch using mapDamage2. Rescaled libraries will not be merged with non-scaled libraries of the same sample for downstream genotyping (this step also often runs into trouble)
10. Bam_trimmed: trimmed ends of mapped reads accordingly (does not use rescaled bam files)
11. Library_merge: merge bam files of same sample but different library type (partial and no UDG together) (if needed)
12. qualimap: report of the trimmed bam files for each sample
13. genotyping: varaint calling of each sample with angsd
14. multiqc: grand summary of all samples

sometimes the damageprofiler and mapdamage_rescaling steps run into problem and the pipeline either exits or get stuck. I decided to remove mapdamage with rescaling, since I am not sure if bam trimming actually works on the rescaled files (I observed bam trimmed files produced before rescaled files). I did rescaling after the pipeline.

 A summary report generated by multiqc at the end of the Eager pipeline can be viewed here for <a href="https://www.dropbox.com/scl/fi/3dajiye78ypi8qd2eplbl/enrichedall_mutliqc.html?rlkey=ocokyrlhmba2p8oii74pb4dif&st=tgjiqpii&dl=0">all target enriched data</a> and
<a href="https://www.dropbox.com/scl/fi/c0gp1mzat82z29e08nj4f/wgsall_multiqc_report.html?rlkey=sdi3xwi0lh5n2kohxh49102gy&st=997v095a&dl=0">all WGS data</a>.

---

## Recalibration
Using trimmed bam as input for mapDamage recalibration (can be considere dthe third step of the eager pipeline)

For enriched data this is done once. For the wgs data, I ran once on the trimmed bam files and another on the additional merged files (which combines udg and non-udg samples together). The latter is for two samples: DVT014 and WMP006. 

[1.0.6_recal_eager.sh](./scripts/1.0.6_recal_eager.sh)
```
#!/bin/bash -l
#SBATCH -J recaleager
#SBATCH --cpus-per-task=8
#SBATCH --time=10-00:00:00
#SBATCH -o /PATH/slurms/slurm-%j-%x.out

# sbatch 1.0.6_recal_eager.sh 00_eager_enrichedall /PATH/genome_HC_allpaths41687_v2.5.fasta /PATH/acrow_eager_list_enrichedall_rescale.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall /PATH/genome_HC_allpaths41687_v2.5.fasta /PATH/acrow_eager_list_wgsall_rescale_2.tsv
# sbatch 1.0.6_recal_eager.sh 00_eager_wgsall /PATH/genome_HC_allpaths41687_v2.5.fasta /PATH/acrow_eager_list_wgsall_rescale_add_2.tsv

echo $(date)
STARTTIME=$(date +%s)

module load charliecloud/0.30

module purge
module load nextflow
conda activate biotools

dat="/PATH"

cd ${dat}
cd $1

nextflow run nf-core/eager -r 2.4.7 -profile conda --input ${3} --fasta ${2} \
-c /PATH/base_modified2.config \
--run_mapdamage_rescaling --skip_preseq --skip_deduplication

#run this again with acrow_eager_list_wgsall_rescale_add_2.tsv to merge the udg and non-udg samples together 

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
```

# Others: Sexing and Mito tree

I also have other scripts to assess sex and species identity with mtDNA. I ran this on the pilot data.

Run [1.0.6_sexid.sh](./scripts/1.0.6_sexid.sh)
```
##sbatch 1.0.6_sexid.sh 00_eager_TE
conda activate biotools
dat="/PATH"

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
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($386/$num<0.6) {print "female_" $386
/$num}else {if ($386/$num>0.9) {print "male_" $386/$num} else {print "unsure_" $3
86/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff16_${i}.txt; done
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($442/$num<0.6) {print "female_" $442
/$num}else {if ($442/$num>0.9) {print "male_" $442/$num} else {print "unsure_" $4
42/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff21_${i}.txt; done
for i in {1..3}; do awk -v num=$i 'NR>1 {if ($563/$num<0.6) {print "female_" $563
/$num}else {if ($563/$num>0.9) {print "male_" $563/$num} else {print "unsure_" $5
63/$num}}}' 01_allchr_cov_transposed.txt > sex_scaff32_${i}.txt; done

paste 01_samples.txt sex_scaff*.txt > 01_sex.txt
```
[1.0.7_mito.sh](./scripts/1.0.7_mito.sh)
```
# sbatch 1.0.7_mito.sh 00_eager_TE
conda activate biotools
module load bcftools
module load seqtk
dat="/PATH"
ref="/PATH/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
mkdir 01_chrM
cd ${1}
cd results/deduplication

for i in *_TE
do
cd $i
# samtools index -b ${i}_rmdup.bam
# samtools view -b ${i}_rmdup.bam chrM > ${i}_chrM.bam
# samtools sort ${i}_chrM.bam > ${i}_chrM_sorted.bam
# samtools index -b ${i}_chrM_sorted.bam -@ 4
# samtools mpileup -f ${ref} -r chrM -Q 20 -q 20 -u ${i}_chrM_sorted.bam |  bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > ${dat}/01_chrM/${i}_chrM.fastq
samtools mpileup -f ${ref} -r chrM -u ${i}_chrM_sorted.bam |  bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > ${dat}/01_chrM/${i}_chrM_relaxed.fastq
# seqtk seq -aQ64 -q20 -n N ${dat}/01_chrM/${i}_chrM.fastq > ${dat}/01_chrM/${i}_chrM.fasta
seqtk seq -aQ64 -q0 -n N ${dat}/01_chrM/${i}_chrM_relaxed.fastq > ${dat}/01_chrM/${i}_chrM_relaxed.fasta
cd ../
done
```
[1.0.7_mito_align.sh](./scripts/1.0.7_mito_align.sh)
```
# sbatch 1.0.7_mito_align.sh 1000
conda activate biotools
module load mafft
module load raxml

dat="/PATH
"
ref="/PATH
/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}
cd 01_chrM

for fname in *_chrM.fasta
do
base=${fname%_chrM.fasta*}
grep -o N ${fname} | wc -l >> wc.tmp 
done

ls *_chrM.fasta > wc2.tmp
paste wc2.tmp wc.tmp > 01_missing_count.txt
rm wc*.tmp

for fname in *_chrM_relaxed.fasta
do
base=${fname%_chrM_relaxed.fasta*}
grep -o N ${fname} | wc -l >> wc_relaxed.tmp
done

ls *_chrM_relaxed.fasta > wc2_relaxed.tmp
paste wc2_relaxed.tmp wc_relaxed.tmp > 01_missing_count_relaxed.txt
rm wc*_relaxed.tmp

#select sequences with less than 8000 Ns (<50% missingness)
for fname in $(awk '{if ($2<8000) print $1}' 01_missing_count.txt)
do
base=${fname%_TE_chrM.fasta*}
awk 'NR >1' ${base}_TE_chrM.fasta | awk -v a=${base} 'BEGIN {print ">"a} {p
rint $0}' >> allmito_selected.fasta
done
mafft --localpair --thread 8 --maxiterate 1000 allmito_selected.fasta > all
mito_selected_linsi_aligned.fasta

#add correctly aligned genbank corvids
cat genbank_realigned.fasta H32_chrM.fasta allmito_selected.fasta > allmito_selected_outgroup.fasta
mafft --localpair --thread 8 --maxiterate 1000 allmito_selected_outgroup.fasta > allmito_selected_outgroup_linsi_aligned.fasta

raxmlHPC -T 8 -s allmito_selected_outgroup_linsi_aligned.fasta -m GTRGAMMA 
-N ${1} -n allmito_selected_outgroup_linsi_aligned_bs${1} -p 12367 -f a -x 16789

#all samples
for fname in $(awk '{if ($2<16000) print $1}' 01_missing_count.txt)
do
base=${fname%_TE_chrM.fasta*}
awk 'NR >1' ${base}_TE_chrM.fasta | awk -v a=${base} 'BEGIN {print ">"a} {pri
nt $0}' >> allmito.fasta
done
mafft --localpair --thread 8 --maxiterate 1000 allmito.fasta > allmito_linsi_
aligned.fasta
cat genbank_realigned.fasta H32_chrM.fasta allmito.fasta > allmito_outgroup.f
asta
mafft --localpair --thread 8 --maxiterate 1000 allmito_outgroup.fasta > allmi
to_outgroup_linsi_aligned.fasta
raxmlHPC -T 8 -s allmito_outgroup_linsi_aligned.fasta -m GTRGAMMA -N ${1} -n 
allmito_outgroup_linsi_aligned_bs${1} -p 12367 -f a -x 16789
```

---
# Variant Calling
For the target enriched data:
51 bones were enriched with two sets of probes produced by myBaits and Twist. I called variants from the combined enriched dataset. The final SNP sets (n=232K) for downstream analysis include ancient and fresh samples.

## Angsd
Call variant by (1) genotype likelihoods, (2) pseudohaploid calling using consensus approach, and (3) absolute genotype across 232K on-target SNP sites. For each approach I processed neutral (n= 194K; 72K outgroup-ascertained transversion) and outlier (n=38K; 1K divergent) SNPs separately.

In the end pseudohaploid variants are used for all analyses, except phenotypic analysis, which requires genotype info and we set a depth of dp>=3

[1.1.1_angsd_combined.sh](./scripts/1.1.1_angsd_combined.sh)
```
#!/bin/bash -l
#SBATCH -J COMangsd
#SBATCH --cpus-per-task=16
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/slurms/slurm-%j-%x.out

# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled neutral probes_193768_neutral.1based
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled outlier probes_38498_outlier.1based
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled neutral_seg72k probes_segAM72k_neutral.1based
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled outlier_1111 knief_snp_1based.pos

#outgroup for admixtool analysis
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled_outgroup neutral probes_193768_neutral.1based #haponly with dumpcount
# sbatch 1.1.1_angsd_combined.sh enrichedallfresh_rescaled_outgroup outlier_1111 knief_snp_1based.pos #haponly dumpcount

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
# awk '{print $1":"$2}' /PATH/01_probes/snp_panel_2/final_snp_panel/neuall_outgroup_ascertained_nooverlap_final_TV_hwe.txt > probes_segAM72k_neutral.1based

cd ${dat}
mkdir 01_angsd_${1}

if [ $(ls ./01_angsd_${1}/bam_*.filelist | wc -l) -eq 1 ]
then
    echo "bam filelist already exist"
else
    echo "generate bam filelist"
    cp bam_${1}.filelist ./01_angsd_${1}/bam_${1}.filelist
fi

cd ${dat}/01_angsd_${1}

missIND=25
minIND=139
#enrichedallfresh_rescaled: 164 ind: 51 ancient + 113 fresh, no dp filtered
#dp3 caused a lot of missingness in anc but no dp cause excussive alt sites

angsd -bam bam_${1}.filelist -dohaplocall 2 -doCounts 1 -rf ${dat}/${3} \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 12 -minMinor 1 -maxMis ${missIND} -out ${1}_${2}_hapconsensus_maxmis_q20 -dumpCounts 2
# then run haploidconversion.ipynb to convert to vcf 

#geno
angsd -bam bam_${1}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} -doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -minInd ${minIND} -remove_bads 1 -geno_minDepth 3 \
-nThreads 12 -rf ${dat}/${3} -out ${1}_${2}_geno_maxmis_q20_dp3

mv ${1}_${2}_geno_maxmis_q20_dp3.tfam ${1}_${2}_geno_maxmis_q20_dp3_original.tfam
awk '{print $2, $1, "0", "0", "0", "-9"}' pop.txt >> ${1}_${2}_geno_maxmis_q20_dp3.tfam
plink --tfile ${1}_${2}_geno_maxmis_q20_dp3 --allow-no-sex --allow-extra-chr -make-bed --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --recode vcf --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --pca --out ${1}_${2}_geno_maxmis_q20_dp3_plink
plink --bfile ${1}_${2}_geno_maxmis_q20_dp3_plink --allow-extra-chr --missing --out ${1}_${2}_geno_maxmis_q20_dp3_plink
bgzip ${1}_${2}_geno_maxmis_q20_dp3_plink.vcf
tabix -p vcf ${1}_${2}_geno_maxmis_q20_dp3_plink.vcf.gz

ENDTIME=$(date +%s)
echo $(date)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
```
---

### PCA of genotype with EMU and SNPrelate
[1.1.3_emu.sh](./scripts/1.1.3_emu.sh)
```
#install emu
#git clone https://github.com/Rosemeis/emu.git
#conda env create -f emu/environment.yml

cd 01_angsd_enrichedallfresh_rescaled

#neutral
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_geno_maxmis_q20_trans_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_geno_maxmis_q20_trans_plink.emu
emu --bfile enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_neutral_seg72k_geno_maxmis_q20_plink.emu

#outlier
emu --bfile enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_dp3_plink.emu
emu --bfile enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_plink --eig 4 --threads 6 --out enrichedallfresh_rescaled_outlier_1111_geno_maxmis_q20_plink.emu
```
Plot using SNPrelate using this [Rscript](./scripts/plot/pca_geno.R).

### PCA of genotype likelihood with PCAngsd
[1.1.2_pcangsd.sh](./scripts/1.1.2_pcangsd.sh)
```
cd ${dat}/01_angsd_${1}
pcangsd --beagle ${1}_${2}.beagle.gz --out ${1}_${2} --threads 8 --snp_weights -pcadapt --tree
```
Plot using this [Rscript](./scripts/plot/pcangsd.R).

### PCA of pseudohaploid with SmartSNP
It is a R package to compute Principal Component Analysis (PCA) on SNP data suitable for ancient and modern DNA (https://christianhuber.github.io/smartsnp). 

The default way of dealing with missing SNP is that if SNP is missing in one or more sample(s) it will be discarded. Alternatively you could choose the imputation method, where missing SNP is imputed with the majority consensus allele type. You could also apply projection method, which project ancient samples in the space of present samples, this method also accomodates for missingness.

This is used for the [main figure 1](./scripts/plot/figure1_map.R). See this Rscripts for only [SmartSNP](./scripts/plot/pca_smartsnp.R).

## Variant calling for the phenotypic analysis
In order to merged with 588 modern genotyped data, we also called genotyped with at least dp=3 across 1111 sites (see [knief_snp_1based.pos](./scripts/knief_snp_1based.pos). 

[1.1.1_angsd_1111snp_AA.sh](./scripts/1.1.1_angsd_1111snp_AA.sh) was run to call genotype from C. moneduloides, which was not included in Knief (2019) as an outgroup.

[1.1.1_angsd_1111snp.sh](./scripts/1.1.1_angsd_1111snp.sh) was run to call genotype on 1111 sites across ancient and fresh samples separately. 
```
conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"

cd ${dat}/02_results
cd ./uli_pca
mkdir 01_angsd_${1}
cd 01_angsd_${1}


#hap_concensus_maxmis_q20
angsd -bam ${dat}/${3}.filelist -dohaplocall 2 -doCounts 1 -rf ../knief_snp_1based.pos \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 \
-nThreads 4 -minMinor 1 -out ${1}_${2}_hapconsensus_q20

#geno
angsd -bam ${dat}/${3}.filelist -GL 2 -doMaf 2 -doMajorMinor 1 -ref ${ref} \
-doGeno 4 -doPost 1 -doGlf 2 -doCounts 1 -doPlink 2 \
-checkBamHeaders 0 -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -geno_minDepth 3 \
-nThreads 4 -rf ../knief_snp_1based.pos -out ${1}_${2}_geno_q20_dp3

#convert to vcf
#first created proper tfam file with proper sample name
mv ${1}_${2}_geno_q20_dp3.tfam ${1}_${2}_geno_q20_dp3_original.tfam
awk '{print $2, $1,"0","0","0","-9"}' pop.txt > ${1}_${2}_geno_q20_dp3.tfam
plink --tfile ${1}_1111_geno_q20_dp3 --allow-no-sex --allow-extra-chr -make-bed --out ${1}_1111_geno_q20_dp3_plink
plink --bfile ${1}_1111_geno_q20_dp3_plink --allow-extra-chr --recode vcf --out ${1}_1111_geno_q20_dp3_plink
plink --bfile ${1}_1111_geno_q20_dp3_plink --allow-extra-chr --pca --out ${1}_1111_geno_q20_dp3_plink

#MODIFY OUTPUT
zcat ${1}_1111_geno_q20_dp3.geno.gz | \
awk '{printf "%s", $1":"$2; for (i=3; i<=NF; i++) {printf "\t%s/%s", substr($i, 1, 1), substr($i, 2, 1);}print "";}' | \
datamash transpose > ${1}_1111_geno_q20_dp3.geno.transposed.txt
```

[Main figure 3](./scripts/plot/figure3_peakheatmap.R) containing genotypes of all 588 modern and 47 (4 removed) ancient samples was plotted using this [Rscript](./scripts/plot/figure3_peakheatmap.R). See below for more details.

---

# Admixtools
First I need to add outgroup moneduloides and jackdaw to the neutral SNPs: [2.2.1_outgroup.sh](./scripts/2.2.1_outgroup.sh)
```
#!/bin/bash -l
#SBATCH -J outgroup
#SBATCH --cpus-per-task=12
#SBATCH --time=14-00:00:00
#SBATCH -o /PATH/05_aDNA/slurms/slurm-%j-%x.out

# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled outlier_geno_maxmis_q20_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh enrichedallfresh_rescaled neutral_geno_maxmis_q20_trans_plink_filtered_norm moneduloides_monedula
# sbatch 2.2.1_outgroup.sh wgs_rescaled all_geno_maxmis_q20_plink_varonly_norepeats moneduloides_monedula

conda activate biotools

dat="/PATH/05_aDNA"
ref="/PATH/05_aDNA/genome_HC_allpaths41687_v2.5.fasta"
out="/PATH/02_outgroups/07_call"

echo $(date)
STARTTIME=$(date +%s)

cd ${dat}/01_angsd_${1}

bcftools query -f '%CHROM\t%POS\n' ${1}_${2}.vcf.gz > ${1}_${2}_positions.txt

#!!!this is used for admixtools!!!
bcftools view --threads 12 -V indels --min-alleles 2 --max-alleles 2 ${1}_${2}_${3}_outgroup.vcf.gz | bgzip -c > ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
tabix -f -p vcf ${1}_${2}_${3}_outgroup_biallele.vcf.gz 
```
Admixtools2 was run on this [Rscript](./scripts/plot/figure2_admixtools_f4_pop_haploid.R) with all neutral pseudohaploid variants. Below shows only the relevant code to generate [main figure 2](./scripts/plot/figure2_admixtools_f4_pop_haploid.R). Refer to [Rscript](./scripts/plot/figure2_admixtools_f4_pop_haploid.R) for all code used to generate admixtools analyses in supplementary.

```
library(magrittr)
library(dplyr)
library(parallel)
library(ape)
library(ggtree, lib.loc="/dss/dsshome1/lxc0E/di67kah/R")

setwd("/PATH/05_aDNA/02_results/admixtools")

neutral = "neutral_all_pop"
seg = "neutral_seg72k_pop"
pop = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k", 
        "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW","MON")
extract_f2(neutral, "neutral_all_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, blgsize = 50000) #allow missing data #160053 SNPs and 27 populations / 57475 for no missing
extract_f2(neutral, "neutral_trans_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, transitions=FALSE, blgsize = 50000) #104164 SNPs and 27 populations / 38422 for no missing
extract_f2(seg, "neutral_seg72k_pop", auto_only = FALSE, pops=pop,overwrite=TRUE,maxmiss=1, blgsize = 50000) #65190 SNPs and 27 populations / 21705 for no missing

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# f4-stats, also known as qpdstat with f4mode=TRUE -----------------------------------------------------------------------------------------------------
# set5: compare EURw to all hoodies across time using irq as pop1
pop1 = c("IRQ0")
pop2 = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k")
pop3 = c("MON")
pop4 = c("JACKDAW","SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k") #jackdaw cannot be in pop2 as it is not closely-related to pop1

rese <- f4(f2_blocks, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
write.table(rese%>%arrange(z), file="./f4est_pseudohaploid/neutral_all_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
rese_trans <- f4(f2_blocks_trans, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
write.table(rese_trans%>%arrange(z), file="./f4est_pseudohaploid/neutral_trans_all_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
rese_seg <- f4(f2_blocks_seg, pop1, pop2, pop3, pop4)%>%
  filter(!pop2 %in% c("EURe100","EURse100"))
write.table(rese_seg%>%arrange(z), file="./f4est_pseudohaploid/neutral_outgroupasc_pop1_irq.txt", quote=FALSE, row.names = FALSE, sep="\t")
#z<3 for all jackdaws and between 20k and SPA

f4st <- rese %>% group_by(pop2) %>% arrange(est, .by_group=TRUE) %>%
  mutate(pop2 = factor(pop2, levels = c("EURs0","RUS0","EURse0","EURsw0","EURn0","EURe0","EURse100","EURe100","EURe1k","EURse1k","NL1k","RUS2k","EURse2k","EURe2k","EURse16k","EURse20k")))
 %>%
  mutate(pop4 = factor(pop4, levels = c("JACKDAW","SPA0","EURw0","EURc0","SPA1k","EURnc1k","EURw2k","EURnw2k","EURc2k","EURc6k"))) 
e1 <- ggplot(f4st%>%filter(pop4!="JACKDAW")) +
  geom_errorbar(aes(y = pop2, x = est, xmin = est - se, xmax = est + se), width = 0.1, alpha=0.5) +
  geom_point(aes(y = pop2, x = est, fill=pop4, shape=pop4), size = 3) +
  scale_fill_manual(values = c("SPA1k"="#7F3B08","EURnc1k"="#FEE0B6", "EURw2k"="#E08214", "EURnw2k"="#B35806", "EURc2k"="#FDB863","EURc6k"="#FDB863")) +
  scale_shape_manual(values = c("JACKDAW"=8,"SPA0"=10,"SPA1k"=21,"EURw0"=3,"EURc0"=4,"EURnc1k"=21, "EURnw2k"=22, "EURw2k"=22, "EURc2k"=22,"EURc6k"=23)) +
  labs(x="f4-estimates (all)", y="Hooded crow") +
  guides(shape=guide_legend("Carrion crow", override.aes = list(shape = c(10,3,4,21,21,22,22,22,23), fill=c("white","white","white","#7F3B08","#FEE0B6","#E08214","#B35806","#FDB863","#FDB
863"))), 
         fill="none") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

# Plot matrix heatmap
me2 <- ggplot(f4st, aes(x = pop4, y = pop2, fill = est, size=est)) +
  geom_point(shape = 21) +
  scale_fill_viridis_c(option = "C",name="f4 estimates") +
  labs(x = "Carrion crow", y = "Hooded crow", size="", fill="") +
  guides(size="none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8))

#plot the simplified graph
tree_text <- "((IRQ,Hooded_crow),(Carrion_crow, Outgroup));"
tree <- read.tree(text = tree_text)
graph2 <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(angle = 0, hjust = 0.5, size = 4)  + coord_flip() + scale_x_reverse()

r1 <- graph2
r2 <- e1 + me2 + plot_layout(widths = c(1,1))
fig3b <- r1 / r2  +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = c("A", "B", "C")) & theme(plot.margin = unit(c(1, 1, 1, 1), "pt"))
fig3b #6x8         

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# qpadm: how target relates to left and right ---------------------------------------------------------------------------------------------------------
# Note: if only 2 outgroup are used, then p-value cannot be generated for rank drop (how well the 2-sources model fit), We rely on pnested to show if 2-sources is significantly better th$
# if >3 outgroups, then p-value is calculated. p_value>0.05 means the model fit well (i cannot reject 2-sources model).
popused = c("SPA1k","SPA0","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k",
            "IRQ0","EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURse1k","EURe1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k","JACKDAW","MON")
f2_blocks = f2_from_precomp("neutral_all_pop", afprod = TRUE, pops=popused) #Discarding 640 block(s)
popa = c("SPA1k","EURw0","EURc0","EURw2k","EURc2k","EURc6k","EURnw2k","EURnc1k")
target=popa
adm1 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm2 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
adm3 <-data.frame(target=popa, SPA0=1, EURse0=0, se=0,z1=0,z2=0,p_rankdrop=0)
left =c("SPA0","EURse0")
right=c("MON","IRQ0","JACKDAW") #no AM for 3 outgroup
#qpadm(f2_blocks, left, right, "SPA1k")

for (i in 1:length(target)) {
  adm1[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm1[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2]
  adm1[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm1[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm1[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm1[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]} #rankdrop$p_nested[1]
write.table(adm1, file="./f4adm_pseudohaploid/neutral_all_target_carrion_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
df_long <- adm1 %>%
  mutate(SPA0 = ifelse(SPA0 < 0, 0, SPA0), EURse0 = ifelse(EURse0 > 1, 1, EURse0)) %>%
  pivot_longer(cols = c(SPA0, EURse0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("SPA1k","EURw0","EURw2k","EURnw2k","EURnc1k","EURc0","EURc2k","EURc6k"))) %>% group_by(target) %>%
  mutate(Mixture = factor(Mixture, levels = c("SPA0","EURse0"))) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)
error <- df_long %>% filter(Mixture == "SPA0") %>% mutate(error_min = 1-ymax - se,error_max = 1-ymax + se)
x1 <- ggplot(df_long, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.05),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) +
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top")

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# qpadm: how target relates to left and right ---------------------------------------------------------------------------------------------------------
popb = c("EURs0","RUS0","EURsw0","EURn0","EURe0","EURse0","EURe100","EURse100","RUS2k","EURe1k","EURse1k","EURe2k","EURse2k","EURse16k","EURse20k","NL1k")
target=popb
adm4 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm5 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
adm6 <-data.frame(target=popb, IRQ0=1, EURc0=0,se=0,z1=0,z2=0,p_rankdrop=0)
left =c("IRQ0","EURc0")
right=c("MON","SPA0","JACKDAW")
#qpadm(f2_blocks_seg, left, right, target)

for (i in 1:length(target)) {
  adm4[i,2] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[1]
  adm4[i,3] <- qpadm(f2_blocks, left, right, target[i])$weights$weight[2]
  adm4[i,4] <- qpadm(f2_blocks, left, right, target[i])$weights$se[1]
  adm4[i,5] <- qpadm(f2_blocks, left, right, target[i])$weights$z[1]
  adm4[i,6] <- qpadm(f2_blocks, left, right, target[i])$weights$z[2]
  adm4[i,7] <- qpadm(f2_blocks, left, right, target[i])$rankdrop$p[1]
  }
write.table(adm4, file="./f4adm_pseudohaploid/neutral_all_target_hooded_3outgrp.txt", quote=FALSE, row.names = FALSE, sep="\t")
df_long4 <- adm4 %>%
  mutate(IRQ0 = ifelse(IRQ0 < 0, 0, IRQ0), EURc0 = ifelse(EURc0 > 1, 1, EURc0)) %>%
  pivot_longer(cols = c(IRQ0, EURc0), names_to = "Mixture", values_to = "Proportion") %>%
  mutate(target = factor(target, levels = c("JACKDAW", "IRQ0","EURs0","EURsw0","EURn0","RUS0","RUS2k","NL1k","EURe0","EURe100","EURe1k","EURe2k","EURse0","EURse100","EURse1k","EURse2k","EURse16k","EURse20k"))) %>% group_by(target) %>%
  arrange(target, Mixture) %>%mutate(ymin = cumsum(lag(Proportion, default = 0)), ymax = ymin + Proportion)
error <- df_long4 %>%filter(Mixture == "EURc0") %>% mutate(error_min = pmax(0, 1 - ymax - se),error_max = pmin(1, 1 - ymax + se))
x2 <- ggplot(df_long4, aes(x = target, y = Proportion, fill = Mixture)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(data = error,aes(ymin = error_min, ymax = error_max),width = 0.4, color = "red") +
  geom_text(data = df_long4,aes(label = paste("p =",sprintf("%.2f", round(p_rankdrop, 2),sep="")), y = 1.01),vjust = 0, size = 3, na.rm = TRUE) +
  scale_fill_manual(values = c("black","grey")) + coord_cartesian(ylim = c(0.00, 1.00))+
  theme_minimal() + labs(x = "Target group (all SNPs)", y = "qpAdm") + theme(legend.position = "top",axis.text.x = element_text(size = 7,angle = 45, vjust = 1, hjust = 1))
```

# Phenotyic analysis on 1111 divergent variants
The variant sites were called using the scripts described above. Genotype heat map across 588 modern and 47 ancient samples were produced using this [Rscript](./scripts/plot/figure3_peakheatmap.R), which is edited from Knief et al., (2019). 

The files required to run this [Rscript](./scripts/plot/figure3_peakheatmap.R) and stored in this [folder](./scripts/02_results/uli_pca/).

First run [3.2_snpeff.sh](./scripts/3.2_snpeff.sh) to annotate the genic nature of each variant 

[3.2_snpeff.sh](./scripts/3.2_snpeff.sh)
```
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
```

## Structural variation
We sequenced ultra long reads using nanopore on 4 selected modern samples, with chr18 enrichment. The scripts for structural variation analysis are in this [folder](./scripts/longreads/). I ran [Sniffles](./scripts/longreads/2.0_sniffles.sh) to look for breakpoint, and [Syri](./scripts/longreads/3.4_syri.sh) to align assembled chr18. Chr18 assembled using [flye assembler](./scripts/longreads/3.0_flye.sh). 

*Note that sample DLoC21 may be referred to as DKoC21 due to naming error from seq company. 

# Selection through time
Ancestry assignment across chromosome based on allele frequencies of reference populations using likelihood approach

We ran this [Rscript](./scripts/plot/figure4_selection.R) to assign each bin of 100 bp region into either Carrion or Hooded crow ancestry. Only delta likelihoods >1 or <-1 are used for confident assignment
