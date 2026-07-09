#---------------------------------------------------------------------------------------------------------------------
conda create -n scikit-env
conda activate scikit-env
conda install -c conda-forge scikit-allel
conda install bioconda::snakemake
cd /projects/psg/people/xcj768/arctic_wolves

#---------------------------------------------------------------------------------------------------------------------
# calc nucleotide diversity in 100kb windows (non-overlapping windows)
cd /projects/psg/people/xcj768/arctic_wolves
conda activate scikit-env
module load samtools/1.21 tabix/1.11 

snakemake -np -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/aw_nt_diversity.smk
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/aw_nt_diversity.arr
#---------------------------------------------------------------------------------------------------------------------
# process pi across chr, calculate estimate of ne
aw_proc_pi.py
#---------------------------------------------------------------------------------------------------------------------
# calc informative windows [proportion of callable sites per window], retain pi estimates in the informative windows, recalculate ne
callable_info_per_window.py
#---------------------------------------------------------------------------------------------------------------------
estimates (for pi & ne) did not change much after filtering by informative windows -
- full dataset [wolf + dog + coyote]: original ne 72,598
filter by windows > 0.7 callable sites: ne 72,343
filter by windows > 0.9 callable sites: ne 69,869

- wolf only dataset, original ne 64,811
filter by windows > 0.7 callable sites: ne 64,498
filter by windows > 0.9 callable sites: ne 62,319
#---------------------------------------------------------------------------------------------------------------------

