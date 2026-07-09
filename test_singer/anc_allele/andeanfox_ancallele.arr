#! /bin/bash
#SBATCH --job-name=aa
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/aa_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/aa_%A-%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G    # memory per cpu-core
#SBATCH --time=8:00:00
#SBATCH --array=1
#SBATCH --qos=normal

#-----------------------------------------------------------------------------------------------------------------------
# extract ancestral allele + annotate vcf
	# could do this for both wolf only & full dataset vcfs - start with wolf only

snakemake -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.smk \
 --use-conda --cores 32 
#-----------------------------------------------------------------------------------------------------------------------