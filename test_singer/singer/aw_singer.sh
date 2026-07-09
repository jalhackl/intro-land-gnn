#---------------------------------------------------------------------------------------------------------------------
#ModuleNotFoundError: No module named 'tskit'
conda create -n tskit-env
conda activate tskit-env
conda install conda-forge::tskit
conda install conda-forge::pandas
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# rewrite to run per chrom [more manageable to control which regions have completed]

cd args/singer/output/test
for i in {1..38}
do
   mkdir -p $i
done

# 1) split vcf per chrom
cd /projects/psg/people/xcj768/arctic_wolves

conda activate tskit-env

module load samtools/1.21 tabix/1.11 
module load  snakemake/9.9.0 
module load gsl/2.5 perl bcftools/1.21  
snakemake  -np -s scripts/downstream/final_dataset/arg/test_singer.smk \
--cores 40 

vi /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/test_singer_smk.arr 
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/test_singer_smk.arr
#---------------------------------------------------------------------------------------------------------------------
# 2) send parallel singer to run per chrom
conda activate tskit-env
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/test_singer1.arr
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# convert parallel singer output to tskit - issues, detailed on github issue https://github.com/popgenmethods/SINGER/issues/39, first advice rerun as regular singer not parallel 
cd /projects/psg/people/xcj768/arctic_wolves
conda activate tskit-env
module load samtools/1.21 tabix/1.11 
module load  snakemake/9.9.0 
module load gsl/2.5 perl bcftools/1.21  
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/converttotskit.arr
# test for chr1, then send for 2-38
# issues here detailed in -  singer_convertotskit_error.sh
#---------------------------------------------------------------------------------------------------------------------
# rerun singer (normal not singer parallel)
# running in screens

# screen 1
conda activate tskit-env
parallel -j 40 '
ID={1}

ENDPOS=$(sed -n "${ID}p" /maps/projects/ilab/people/xcj768/arctic_wolf/args/canFam31_autosomes_lengths.bed | awk "{print \$3+1}")

echo "${ID} -> ${ENDPOS}"

/maps/projects/ilab/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/singer_master \
  -Ne 72343 \
  -m 4.96e-9 \
  -recomb_map /maps/projects/ilab/people/xcj768/arctic_wolf/R_maps/${ID}.average_canFam3.1.txt \
  -polar 0.99 \
  -vcf /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/input/aw_n82_${ID} \
  -start 0 \
  -end ${ENDPOS} \
  -output /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/output/${ID}/aw_n82_${ID}' ::: {1..2}

# screen 2
conda activate tskit-env
parallel -j 80 '
ID={1}

ENDPOS=$(sed -n "${ID}p" /maps/projects/ilab/people/xcj768/arctic_wolf/args/canFam31_autosomes_lengths.bed | awk "{print \$3+1}")

echo "${ID} -> ${ENDPOS}"

/maps/projects/ilab/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/singer_master \
  -Ne 72343 \
  -m 4.96e-9 \
  -recomb_map /maps/projects/ilab/people/xcj768/arctic_wolf/R_maps/${ID}.average_canFam3.1.txt \
  -polar 0.99 \
  -vcf /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/input/aw_n82_${ID} \
  -start 0 \
  -end ${ENDPOS} \
  -output /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/output/${ID}/aw_n82_${ID}' ::: {3..38}

screen -r test3_singer

# convert to trees
cd /maps/projects/ilab/people/xcj768/arctic_wolf/
cd /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/output/38


ID=38
input="aw_n82_${ID}"
outp="trees/aw_n82_${ID}"
python /maps/projects/ilab/people/xcj768/arctic_wolf/scripts/args/convert_to_tskit_nspope.py \
  -input ${input} \
  -output ${outp} \
  -start 0 \
  -end 99

 screen -r test_singer

