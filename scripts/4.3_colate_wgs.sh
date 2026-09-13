#!/bin/bash -l
#SBATCH -J colate
#SBATCH --mail-user=gwee@biologie.uni-muenchen.de
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --exclude=hlegr1409n06,hlegr1409n07
#SBATCH -o PATH/05_aDNA/slurms/slurm-%j-%x.out

#sbatch 4.3_colate_wgs.sh BRW001 all_chr 1873 E02 hooded  all_chr
#sbatch 4.3_colate_wgs.sh DVT014 all_chr 16047 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh DVT017 all_chr 16007 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh KCZ003 all_chr 428 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh KCZ012 all_chr 724 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh NCP001 all_chr 1545 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh VKP001 all_chr 1176 E02 hooded all_chr
#sbatch 4.3_colate_wgs.sh SLN001 all_chr 1178 E02 carrion all_chr
#sbatch 4.3_colate_wgs.sh TDN002 all_chr 910 E02 carrion all_chr
#sbatch 4.3_colate_wgs.sh WMP006 all_chr 1926 E02 carrion all_chr
#sbatch 4.3_colate_wgs.sh E02 peak_chr18 0 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh BRW001 peak_chr18 1873 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh DVT014 peak_chr18 16047 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh DVT017 peak_chr18 16007 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh KCZ003 peak_chr18 428 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh KCZ012 peak_chr18 724 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh NCP001 peak_chr18 1545 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh VKP001 peak_chr18 1176 E02 hooded chr18
#sbatch 4.3_colate_wgs.sh SLN001 peak_chr18 1178 E02 carrion chr18
#sbatch 4.3_colate_wgs.sh TDN002 peak_chr18 910 E02 carrion chr18
#sbatch 4.3_colate_wgs.sh WMP006 peak_chr18 1926 E02 carrion chr18

echo $(date)
STARTTIME=$(date +%s)

conda activate biotools
dat="PATH/05_aDNA"
ref="PATH/06_longreads/ASM73873v6"

cd ${dat}
cd ./02_results/colate
mkdir ./wgs/${2}
#git clone https://github.com/leospeidel/Colate.git
#chmod +x binaries/v0.1.5_x86_64_dynamic/bin/Colate

#generate colate input files
PATH/05_aDNA/Colate/binaries/v0.1.5_x86_64_dynamic/bin/Colate \
	--mode make_tmp \
	--mut ./data/mut-ages/relate \
	--target_bam ${dat}/00_eager_wgsall_02/results/mapped_ref5.7/${1}_ref5.7_${6}.bam \
	--ref_genome ${ref} \
	    --chr ${2}.txt \
		--filters 20,20,20 \
	-o ./wgs/${2}/out_${1}_${2}

#compute coalescence rates THIS WILL BE REGENERATED AGAIN WITH THE PYTHON SCRIPT (run_colate_pairwise.py)
bins="3,7,0.2" #epochs in log10 years (format: start,end,stepsize)
/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2020__ancientDNA/05_aDNA/Colate/binaries/v0.1.5_x86_64_dynamic/bin/Colate \
	--mode mut \
	--mut ./data/mut-ages/relate \
	--target_tmp ./wgs/${2}/out_${1}_${2}.colate.in \
	--reference_tmp  ./wgs/${2}/out_${4}_${2}.colate.in \
        --bins ${bins} \
	--chr ${2}.txt \
	--num_bootstraps 20 \
	--target_age ${3} \
	--reference_age 0 \
	--years_per_gen 5.79 \
	-o ./wgs/${2}/coalescence_${5}_${4}_${1}_${2}


ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
