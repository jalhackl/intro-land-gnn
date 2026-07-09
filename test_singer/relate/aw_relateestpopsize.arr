#! /bin/bash
#SBATCH --job-name=relateestpopsize
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/runrelateestpopsize_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/runrelateestpopsize_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=120:00:00
#SBATCH --array=1
#SBATCH --qos=normal
#-----------------------------------------------------------------------------------------------------------------------
# estimate popsize
chrom=$SLURM_ARRAY_TASK_ID
path_to_relate="/projects/psg/people/xcj768/resources/relate_v1.2.4_x86_64_dynamic"
mu=4.96e-9
gen=3
pref="aw_n81"
input_poplabels="/projects/psg/people/xcj768/arctic_wolves/args/relate/input/aw.poplabels"
#prefix=${pref}"_"${chrom}


${path_to_relate}"/scripts/EstimatePopulationSize/EstimatePopulationSize.sh" \
    -i ${pref} \
    -m ${mu} \
    --years_per_gen ${gen} \
    --poplabels ${input_poplabels} \
    --threads 8 \
    --first_chr 1 \
    --last_chr 38 \
    -o ${pref}"_popsize"
#-----------------------------------------------------------------------------------------------------------------------
