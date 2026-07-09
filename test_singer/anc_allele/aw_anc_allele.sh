#-----------------------------------------------------------------------------------------------------------------------
# get ancestral allele from outgroup [only homozygous calls]
# annotate AA in the imputed + filtered bcf

cd /projects/psg/people/xcj768/arctic_wolves
module load gsl/2.5 perl bcftools/1.21  
module load openjdk/20.0.0 gcc/13.2.0 R/4.4.2
module load samtools/1.21 tabix/1.11 

snakemake -np -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/andeanfox_ancallele.smk

sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.arr
#-----------------------------------------------------------------------------------------------------------------------
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n'  anc_allele/aw_infosc0.8_maf0.01.autosomes.wolf_n56_aa.vcf.gz | head 
chr1	1763	T	G	AA%3DT
chr1	2584	T	A	AA%3DT
chr1	2622	C	G	AA%3DC

#%3D is URL encoding for the character =
# ancestral allele has annotated, but instead of = is %3D

#-----------------------------------------------------------------------------------------------------------------------
# then will need to polarise - flip where alt is the anc state
# flipping alleles in r too heavy -> python
#-----------------------------------------------------------------------------------------------------------------------
conda create -n cyvcf2-env
conda activate cyvcf2-env
conda install bioconda::cyvcf2
conda install bioconda::snakemake
#-----------------------------------------------------------------------------------------------------------------------
cd /projects/psg/people/xcj768/arctic_wolves
conda activate cyvcf2-env
conda install bioconda::snakemake

module load gsl/2.5 perl bcftools/1.21  
module load openjdk/20.0.0 gcc/13.2.0 R/4.4.2
module load samtools/1.21 tabix/1.11 

snakemake -np -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.smk

less /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.arr
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.arr
#-----------------------------------------------------------------------------------------------------------------------
# index the vcfs
cd /projects/psg/people/xcj768/arctic_wolves
module load gsl/2.5 perl bcftools/1.21  

bcftools index -f anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.wolf_n56_aa_polarised.vcf.gz
bcftools index -f anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.n82.nooutgroup_aa_polarised.vcf.gz

#-----------------------------------------------------------------------------------------------------------------------
# check allele flipping has been successful
check_allele_flips.R
# has not worked => need to rewrite & resend

# resubmitted 
cd /projects/psg/people/xcj768/arctic_wolves
conda activate cyvcf2-env

module load gsl/2.5 perl bcftools/1.21  

module load openjdk/20.0.0 gcc/13.2.0 R/4.4.2
module load samtools/1.21 tabix/1.11 

snakemake -np -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.smk
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/andeanfox_ancallele.arr
#-----------------------------------------------------------------------------------------------------------------------
# check allele flips - has now worked
check_allele_flips.R
#-----------------------------------------------------------------------------------------------------------------------


