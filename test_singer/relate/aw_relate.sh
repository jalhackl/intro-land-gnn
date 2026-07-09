# transfer relate binary to cluster

#-----------------------------------------------------------------------------------------------------------------------
# make --poplabels Optional: File containing population labels of samples.
mk_relate_poplabels.R

#-----------------------------------------------------------------------------------------------------------------------
# mk input files for relate

cd /projects/psg/people/xcj768/arctic_wolves
conda activate glimpse-env
module load gsl/2.5 perl bcftools/1.21  
module load samtools/1.21 tabix/1.11 

snakemake -np -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
 --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml

sbatch  /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/mk_relate_input.arr 
#-----------------------------------------------------------------------------------------------------------------------
# mk mask fasta

# make inverse of good regions - ie regions to mask out
awk '{print $1"\t"$3}' site_filters/depth/canFam31_autosomes_lengths.bed > args/relate/input/mk_mask/canFam31_autosomes_lengths.txt
cd /projects/psg/people/xcj768/arctic_wolves
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/mk_relatemask.arr
# checked makes sense
# fmt for relate - Loci passing the mask are denoted by P, loci not passing the mask are denoted by N.
sed -r 's/[ACGT]/P/g' args/relate/input/mk_mask/goodsites_mask.fa > args/relate/input/mk_mask/goodsites_mask_relatefmt.fa 
#-----------------------------------------------------------------------------------------------------------------------
# mk input files for relate, include mask -> segmentation fault here, need to troubleshoot for now generate without mask
cd /projects/psg/people/xcj768/arctic_wolves
conda activate glimpse-env
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/mkrelateinp_plusmask.arr
#-----------------------------------------------------------------------------------------------------------------------
cd /projects/psg/people/xcj768/arctic_wolves
conda activate glimpse-env
module load gsl/2.5 perl bcftools/1.21  
module load samtools/1.21 tabix/1.11 

snakemake -np -s scripts/downstream/final_dataset/relate/mk_relate_input.smk \
 --configfile scripts/downstream/final_dataset/relate/mk_relate_input.yaml
vi /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/mk_relate_input.arr 
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/mk_relate_input.arr 
# input files generated
#-----------------------------------------------------------------------------------------------------------------------
# mk soft links to recomb maps per chrom
mkdir -p args/relate/input/R_maps
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/input/R_maps
for ID in {1..38};
do
ln -s /projects/psg/people/qvw641/ResQ_wolf/HapNe/VCFs/${ID}.average_canFam3.1.txt 
done
#-----------------------------------------------------------------------------------------------------------------------
# prep to run relate mode all
-N,--effective_size
Effective population size of haplotypes. (NOT of individuals! To get the population size of haplotypes, multiply the effective population size of individuals by 2)
> 72343*2
[1] 144686

# mk dictionary for yaml (paths to recombination maps)
for ID in {1..38};
do
echo '  "'${ID}'": /projects/psg/people/qvw641/ResQ_wolf/HapNe/VCFs/'${ID}'.average_canFam3.1.txt' 
done

#-----------------------------------------------------------------------------------------------------------------------
# run relate mode all [needs to be run in output dir] -> this had mem issues (died with <Signals.SIGKILL: 9>.) go to next section
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/mode_all 
conda activate glimpse-env
module load gsl/2.5 perl bcftools/1.21  
module load samtools/1.21 tabix/1.11 

snakemake -np -s /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/run_relate_rule_all.smk \
 --configfile /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/run_relate_all.yaml

sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/run_relate_all.arr
# Submitted batch job 44516498  # Sat May 16 12:13:12 CEST 2026
#-----------------------------------------------------------------------------------------------------------------------
# run relate parallel
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output
conda activate glimpse-env
module load gsl/2.5 perl bcftools/1.21
module load samtools/1.21 tabix/1.11
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/aw_relateparallel.arr
#-----------------------------------------------------------------------------------------------------------------------
# estimate population sizes, expects input files to be named pref_chr1, etc
cd  /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
for ID in {1..37}
do
ln -s ../aw_n81_${ID}.mut aw_n81_chr${ID}.mut
ln -s ../aw_n81_${ID}.anc aw_n81_chr${ID}.anc
done

conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/aw_relateestpopsize.arr 
#-----------------------------------------------------------------------------------------------------------------------
# install helper lib to convert relate to tskit format
install_relatelib.sh

# convert relate output (of all pops estimatepopulationsize) [format anc, mut] to tskit format [trees]
module load gcc/10.2.0 
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize/
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/relatetotskit.arr 
#-----------------------------------------------------------------------------------------------------------------------
cd /projects/psg/people/xcj768/arctic_wolves/
conda env list
conda activate tskit-env
conda list
conda install conda-forge::tskit
#The following packages will be UPDATED:
#  ca-certificates                       2026.1.4-hbd8a1cb_0 --> 2026.5.20-hbd8a1cb_0 
#  openssl                                  3.6.0-h26f9b46_0 --> 3.6.2-h35e630c_0 
#  tskit                               1.0.0-py314hc02f841_0 --> 1.0.3-py314hc02f841_0 
conda install conda-forge::tszip
#-----------------------------------------------------------------------------------------------------------------------
# add metadata to relate trees, following https://github.com/leospeidel/relate_lib/issues/7
cd  /projects/psg/people/xcj768/arctic_wolves/
conda activate tskit-env
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/addmetadata_relatetrees.arr
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# selection in greenland/ellesmere wolves only - first estimate the population size history using this module, with setting --threshold 0 - go to next section
conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/grell_relateestpopsize.arr 
# error here - 
(base) [xcj768@mjolnircomp02fl arctic_wolves]$ tail  /projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/rungrellrelateestpopsize_44714368-1.err
Combining .anc/.mut files into one file...
CPU Time spent: 22.116243s; Max Memory usage: 4.468Mb.
---------------------------------------------------------

---------------------------------------------------------
Finalizing coalescence rate...
22 162
Error: number of haplotypes in anc/mut does not match number of samples in .poplabels file
You can just rerun this step using:
PATH_TO_RELATE/bin/RelateCoalescentRate --mode FinalizePopulationSize -o example --poplabels example.poplabels

#-----------------------------------------------------------------------------------------------------------------------
# try to use pop of itnerest flag - this works [greenland, ellesmere wolves threshold 0], then next step detect selection
conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/grell_relateestpopsize1.arr 
#-----------------------------------------------------------------------------------------------------------------------
# detect selection in greenland+ellesmere wolves

conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/grell_relateseln.arr

#-----------------------------------------------------------------------------------------------------------------------
# all arctic wolves threhsold 0 [step pre detect selection]
conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/arctics_relateestpopsize1.arr 

# detect selection - all arctic wolves
conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize/all_arctics
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/arctics_relateseln.arr 

#-----------------------------------------------------------------------------------------------------------------------
# run estpopsize per popn

array=(wEllesmere wGreenland wNorthBaffin wSouthBaffin wVictoriaIsland wBanksIsland wCentral wAlaska wAlberta dAlaska dGreenland)
for value in "${array[@]}"
do
   echo $value >> /projects/psg/people/xcj768/arctic_wolves/args/relate/input/aw.popsize.popns.txt
done

array=(wEllesmere wGreenland wNorthBaffin wSouthBaffin wVictoriaIsland wBanksIsland wCentral wAlaska wAlberta dAlaska dGreenland)
for value in "${array[@]}"
do
   mkdir -p "per_popn/"$value
done


# run estpopsize per popn
conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/perpop_relatestpopsize.arr

#-----------------------------------------------------------------------------------------------------------------------

# greenland+ellesmere wolves sample branch lengths - this needs to be done for individual snps [to be prioritised] [this step has not been run]

conda activate r_env
cd /projects/psg/people/xcj768/arctic_wolves/args/relate/output/popsize
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/relate/grell_relatesamplebranchlengths.arr
#-----------------------------------------------------------------------------------------------------------------------

