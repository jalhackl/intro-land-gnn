#! /bin/bash
#SBATCH --job-name=ntdiversity
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/ntdiversity_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/ntdiversity_%A-%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G    # memory per cpu-core
#SBATCH --time=8:00:00
#SBATCH --array=1
#SBATCH --qos=normal

#-----------------------------------------------------------------------------------------------------------------------
# calc nucleotide diversity in 100kb windows
snakemake -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/aw_nt_diversity.smk \
 --use-conda --cores 32 
#-----------------------------------------------------------------------------------------------------------------------