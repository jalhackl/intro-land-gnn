#! /bin/bash
#SBATCH --job-name=relateinput
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/mkrelateinput_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/mkrelateinput_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH --array=1
#SBATCH --qos=normal
#-----------------------------------------------------------------------------------------------------------------------
# step done
#bcftools consensus  -f gts/filter_ref_fasta/canFam31_autosomes.fasta \
# -s AndeanFox impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.n82.bcf > args/relate/input/cafam31_anc.fa

# MissingOutputException in rule mk_relate_input - check this
#snakemake -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
# --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml \
#--use-conda --cores 32 --unlock

#snakemake -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
# --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml \
#--use-conda --cores 32 --rerun-incomplete
#-----------------------------------------------------------------------------------------------------------------------
# testing -
#/projects/psg/people/xcj768/resources/relate_v1.2.4_x86_64_dynamic/scripts/PrepareInputFiles/PrepareInputFiles.sh \
# --haps args/relate/input/haps/aw_n81_8.haps \
# --sample args/relate/input/haps/aw_n81_8.sample \
# --ancestor args/relate/input/cafam31_anc.fa  \
# --poplabels args/relate/input/aw.poplabels \
# -o args/relate/input/final/aw_n81_8
#-----------------------------------------------------------------------------------------------------------------------
snakemake -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
 --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml \
--use-conda --cores 32 --unlock

snakemake -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
 --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml \
--use-conda --cores 32 --rerun-incomplete

