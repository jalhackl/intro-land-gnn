#---------------------------------------------------------------------------------------------------------------------
import argparse
import allel
import pickle
#-------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="input params from smk")
parser.add_argument("--input_vcf", type=str, help="input vcf")
parser.add_argument("--chrom", type=str, help="chrom")
parser.add_argument("--pref", type=str, help="pref")
#-------------------------------------------------------------------------------------------------
args = parser.parse_args()
input_vcf = args.input_vcf
chrom = args.chrom
pref = args.pref
#---------------------------------------------------------------------------------------------------------------------
# read in per chrom
callset = allel.read_vcf(input_vcf, region=chrom, fields=['CHROM','POS','GT'])
pos = callset['variants/POS']
gt = allel.GenotypeArray(callset['calldata/GT'])
ac = gt.count_alleles()
# calc nucleotide diversity in windows, size 100kb
pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=100000)
with open(f"genetic_diversity/nucl_diversity/{pref}_{chrom}_pi", "wb") as f:
        pickle.dump((pi, windows, n_bases, counts), f)
#---------------------------------------------------------------------------------------------------------------------