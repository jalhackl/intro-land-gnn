import pandas as pd
import argparse
from cyvcf2 import VCF, Writer
#-------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="input params from smk")
parser.add_argument("--input_vcf", type=str, help="input vcf")
parser.add_argument("--output_vcf", type=str, help="output vcf")
#-------------------------------------------------------------------------------------------------
args = parser.parse_args()
input_vcf = args.input_vcf
output_vcf = args.output_vcf
#-------------------------------------------------------------------------------------------------
anc = pd.read_csv('anc_allele/andeanfox_anc.gz',sep='\t', header=None)
anc.columns=['chr','pos','ref','alt','aa']
#            chr       pos ref alt    aa
#0          chr1      1763   T   G  AA=T
#1          chr1      2584   T   A  AA=T
#2          chr1      2622   C   G  AA=C
#3          chr1      2856   G   A  AA=G
#4          chr1      3147   G   A  AA=G
#...         ...       ...  ..  ..   ...
#15578635  chr38  23914372   A   G  AA=G
#15578636  chr38  23914426   C   A  AA=A
#15578637  chr38  23914464   C   A  AA=C
#15578638  chr38  23914482   C   A  AA=.
#15578639  chr38  23914497   G   A  AA=.
#[15578640 rows x 5 columns]
# split anc allele field
anc[['aa','al']] = anc['aa'].str.split('AA=',expand=True)
#>>> anc
#            chr       pos ref alt aa al
#0          chr1      1763   T   G     T
#1          chr1      2584   T   A     T
#2          chr1      2622   C   G     C
#-------------------------------------------------------------------------------------------------
anc_dict = anc.set_index(["chr", "pos"])["al"].to_dict()
#-------------------------------------------------------------------------------------------------
# check 1 site that needs to be flipped
#7  chr1 3500  A  G AA=G        G
#from cyvcf2 import VCF, Writer
#input_vcf="anc_allele/aw_infosc0.8_maf0.01.autosomes.wolf_n56_aa.vcf.gz"
#vcf = VCF(input_vcf)
#for record in vcf("chr1:3500-3501"):
#    print("REF:", record.REF, type(record.REF))
#    print("ALT:", record.ALT, type(record.ALT))
#REF: A <class 'str'>
#ALT: ['G'] <class 'list'>
#old_ref = record.REF
#old_alt = record.ALT[0]
#record.REF = old_alt
#record.ALT = [old_ref] 
#>>> old_ref
#'A'
#>>> old_alt
#'G'
#>>> record.REF
#'G'
#>>> record.ALT
#['A']
#new_gts = []
#for a, b, phased in record.genotypes:
#    new_gts.append([
#        1 if a == 0 else 0 if a == 1 else a,
#        1 if b == 0 else 0 if b == 1 else b,
#        phased
#])
#>>> record.genotypes
#[[1, 0, True], [0, 0, True], [0, 0, True], [1, 0, True], [1, 0, True], [0, 0, True], [0, 0, True], [0, 1, True], [1, 0, True], [1, 0, True], [0, 0, True], [1, 0, True], [0, 0, True], [0, 0, True], [1, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [1, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [0, 0, True], [0, 1, True], [1, 0, True], [1, 0, True], [0, 1, True], [1, 0, True], [1, 0, True], [0, 1, True], [0, 1, True], [0, 0, True], [0, 0, True], [1, 0, True], [0, 0, True], [1, 0, True], [1, 0, True], [1, 0, True], [0, 0, True], [0, 0, True], [1, 0, True], [0, 0, True], [1, 0, True], [1, 1, True], [1, 0, True], [1, 0, True], [1, 0, True], [0, 1, True], [1, 0, True], [1, 0, True], [1, 0, True], [1, 0, True], [1, 0, True]]
#>>> new_gts
#[[0, 1, True], [1, 1, True], [1, 1, True], [0, 1, True], [0, 1, True], [1, 1, True], [1, 1, True], [1, 0, True], [0, 1, True], [0, 1, True], [1, 1, True], [0, 1, True], [1, 1, True], [1, 1, True], [0, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [0, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [1, 1, True], [1, 0, True], [0, 1, True], [0, 1, True], [1, 0, True], [0, 1, True], [0, 1, True], [1, 0, True], [1, 0, True], [1, 1, True], [1, 1, True], [0, 1, True], [1, 1, True], [0, 1, True], [0, 1, True], [0, 1, True], [1, 1, True], [1, 1, True], [0, 1, True], [1, 1, True], [0, 1, True], [0, 0, True], [0, 1, True], [0, 1, True], [0, 1, True], [1, 0, True], [0, 1, True], [0, 1, True], [0, 1, True], [0, 1, True], [0, 1, True]]
# this now looks correct
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
# read in vcf
vcf = VCF(input_vcf)
w = Writer(output_vcf, vcf)
#-------------------------------------------------------------------------------------------------
# iterate over vcf
    # flip ref/alt where alt is ancestral
        # flip alleles
        # flip gts
for record in vcf:
    key = (record.CHROM, record.POS)
    ancestral = anc_dict.get(key)
    if (ancestral not in (None, ".") and record.ALT and len(record.ALT) == 1 and record.ALT[0] == ancestral):
        old_ref = record.REF
        old_alt = record.ALT[0]
        record.REF = old_alt
        record.ALT = [old_ref]
        new_gts = []
        for a, b, phased in record.genotypes:
            new_gts.append([
                1 if a == 0 else 0 if a == 1 else a,
                1 if b == 0 else 0 if b == 1 else b,
                phased])
        record.genotypes = new_gts
    w.write_record(record)
w.close()
#-------------------------------------------------------------------------------------------------