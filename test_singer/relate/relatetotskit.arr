#! /bin/bash
#SBATCH --job-name=relatetotskit
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/converttotskit_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/converttotskit_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --array=1-37
#SBATCH --qos=normal
#-----------------------------------------------------------------------------------------------------------------------
# convert from Relate to tskit format

# first one chrom (38), then 1-37

chrom=$SLURM_ARRAY_TASK_ID
path_to_relate_lib="/projects/psg/people/xcj768/resources/relate_lib"

${path_to_relate_lib}/bin/Convert --mode ConvertToTreeSequence \
              --compress \
              --anc aw_n81_popsize_chr${chrom}.anc.gz \
              --mut aw_n81_popsize_chr${chrom}.mut.gz \
              -o tskit/aw_n81_popsize_chr${chrom}
#-----------------------------------------------------------------------------------------------------------------------
# compresses these Relate-converted tree sequences by assigning the same age to nodes with identical descendant sets across adjacent trees.