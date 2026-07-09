#! /bin/bash
#SBATCH --job-name=relateparallel
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/runrelateparallel_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/runrelateparallel_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=18:00:00
#SBATCH --array=1-38
#SBATCH --qos=normal
#-----------------------------------------------------------------------------------------------------------------------
# parallel relate
chrom=$SLURM_ARRAY_TASK_ID
path_to_relate="/projects/psg/people/xcj768/resources/relate_v1.2.4_x86_64_dynamic"
mu=4.96e-9
haplotype_n=144686
pref="aw_n81"
input_haps="/projects/psg/people/xcj768/arctic_wolves/args/relate/input/final/"${pref}"_"${chrom}".haps.gz"
input_sample="/projects/psg/people/xcj768/arctic_wolves/args/relate/input/final/"${pref}"_"${chrom}".sample.gz"
input_annot="/projects/psg/people/xcj768/arctic_wolves/args/relate/input/final/"${pref}"_"${chrom}".annot"
input_rmap="/projects/psg/people/qvw641/ResQ_wolf/HapNe/VCFs/"${chrom}".average_canFam3.1.txt"
prefix=${pref}"_"${chrom}


${path_to_relate}"/scripts/RelateParallel/RelateParallel.sh" \
-m ${mu} \
-N ${haplotype_n} \
--haps ${input_haps} \
--sample ${input_sample} \
--map ${input_rmap} \
--annot ${input_annot} \
-o ${prefix} \
--threads 8 
#-----------------------------------------------------------------------------------------------------------------------

