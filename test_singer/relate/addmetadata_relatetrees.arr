#! /bin/bash
#SBATCH --job-name=relateaddmeta
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/addrelatemetad_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/addrelatemetad_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --array=1-37
#SBATCH --qos=normal
#-----------------------------------------------------------------------------------------------------------------------
# addmetadata per chrom to relate trees
chrom=$SLURM_ARRAY_TASK_ID

script_dir="scripts/downstream/final_dataset/relate"
inputtree="args/relate/output/popsize/tskit/aw_n81_popsize_chr"${chrom}".trees"
poplabels="args/relate/input/aw.poplabels"
out_prefix="args/relate/output/popsize/tskit/aw_n81_popsize_plmetadata_chr"${chrom}

python ${script_dir}/addmetadata_relatetrees.py \
        --ts_file ${inputtree} \
        --poplabels_file ${poplabels} \
        --out_prefix ${out_prefix}

#-----------------------------------------------------------------------------------------------------------------------