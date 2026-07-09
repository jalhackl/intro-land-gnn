#! /bin/bash
#SBATCH --job-name=parallelsinger
#SBATCH --output=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/parallelsinger_%A-%a.log
#SBATCH --error=/projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/parallelsinger_%A-%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --array=1-38
#SBATCH --qos=normal

#-----------------------------------------------------------------------------------------------------------------------
ID=$SLURM_ARRAY_TASK_ID
#-----------------------------------------------------------------------------------------------------------------------
# run parallel singer
module load parallel/20241222 

input="/projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_${ID}"
outp="/projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/${ID}/aw_n82_${ID}"

/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/parallel_singer \
-Ne 72343 \
-m 4.96e-9 \
-vcf ${input} \
-polar 0.99 \
-output ${outp} \
-n 100 \
-thin 20

log="/projects/psg/people/xcj768/arctic_wolves/args/singer/logs/${ID}.log"
touch ${log}
echo 'singer finished for chrom ${ID}' >> ${log}
#---------------------------------------------------------------------------------------------------------------------
#options:
#  -h, --help            show this help message and exit
#  -Ne NE                Effective population size. Default: 1e4.
#  -m M                  Mutation rate.
#  -ratio RATIO          Recombination to mutation ratio. Default: 1.
#  -L L                  Block length. Default: 1e6.
#  -vcf VCF              VCF file prefix (without .vcf or .vcf.gz extension).
#  -output OUTPUT        Output file prefix.
#  -n N                  Number of MCMC samples.
#  -thin THIN            Thinning interval length.
#  -polar POLAR          Site flip probability. Default: 0.5.
#  -freq FREQ            Convert to tskit every {freq} samples. Default: 1.
#  -num_cores NUM_CORES  Number of cores. Default: 20.
#---------------------------------------------------------------------------------------------------------------------
