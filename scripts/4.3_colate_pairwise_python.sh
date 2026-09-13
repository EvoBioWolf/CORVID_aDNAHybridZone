#!/bin/bash -l
#SBATCH -J colate_pair
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH -o PATH/05_aDNA/slurms/slurm-%j-%x.out

#sbatch 4.3_colate_pairwise_python.sh all_chr wgs
#sbatch 4.3_colate_pairwise_python.sh peak_chr18 wgs

#sbatch 4.3_colate_pairwise_python.sh all_chr enriched
#sbatch 4.3_colate_pairwise_python.sh peak_chr18 enriched

conda activate renv

dat="PATH/05_aDNA"

python_script="${dat}/02_results/colate/heatmap/run_colate_pairwise.py"
TYP=$2
CHR=$1
N_WORKERS=8

cd ${dat}
python3 "${python_script}" "${CHR}" "${TYP}" "${N_WORKERS}"

ENDTIME=$(date +%s)
echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task"
