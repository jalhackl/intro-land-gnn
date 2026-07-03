import os
import pandas as pd
#import pybedtools
import numpy as np

def list_subdirectories(directory):
    subdirectories = []
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            subdirectories.append(item)
    return subdirectories



def ms2vcf(ms_file, vcf_file, nsamp, seq_len, ploidy=2, ind_prefix="ms_"):
    """
    Description:
        Converts ms output files into the VCF format.

    Arguments:
        ms_file str: Name of the ms file (input).
        vcf_file str: Name of the VCF file (output).
        nsamp int: Number of haploid genomes.
        seq_len int: Sequence length.
        ploidy int: Ploidy of each individual.
    """
    data = []
    i = -1
    header = "##fileformat=VCFv4.2\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([ind_prefix + str(i) for i in range(int(nsamp/ploidy))])

    with open(ms_file, 'r') as f:
        #f.readline()
        #f.readline()
        next(f)
        next(f)
        for l in f:
        #for l in f.readlines():
            if l.startswith('//'):
                i += 1
                data.append({})
                data[i]['pos'] = []
                data[i]['geno'] = []
            elif l.startswith('positions'):
                data[i]['pos'] = l.rstrip().split(" ")[1:]
            elif l.startswith('0') or l.startswith('1'):
                data[i]['geno'].append(l.rstrip())

    shift = 0
    with open(vcf_file, 'w') as o:
        o.write(header+"\n")
        for i in range(len(data)):
            for j in range(len(data[i]['pos'])):
                pos = int(seq_len * float(data[i]['pos'][j])) + shift
                genotypes = "".join([data[i]['geno'][k][j] for k in range(len(data[i]['geno']))])
                genotypes = "\t".join([a+'|'+b for a,b in zip(genotypes[0::ploidy],genotypes[1::ploidy])])
                o.write(f"1\t{pos}\t.\tA\tT\t100\tPASS\t.\tGT\t{genotypes}\n")
            shift += seq_len



def create_ind_lists(nsamp, ss_ref_list, ss_tgt_list, ploidy=2, ind_prefix="ms_"):
    nr_of_inds = (int(nsamp/ploidy))
    ntgt = 1
    nref = nr_of_inds - ntgt
    with open(ss_tgt_list, 'w') as f:
        for i in range(ntgt):
            f.write(f'{ind_prefix}{i}\n')
    with open(ss_ref_list, 'w') as f:
        for i in range(i+1, nref + ntgt):
            f.write(f'{ind_prefix}{i}\n')


# fn  process_sstar_1src_output moved to workflow/scripts/sstar_sprime_out_to_bed.py
def process_sstar_1src_output(in_file, out_file):
    """
    Description:
        Helper function for converting output from sstar to BED format files.

    Arguments:
        in_file str: Name of the input file.
        out_file str: Name of the output file.
    """
    df = pd.read_csv(in_file, sep="\t")
    df = df[df['significant'] == True]
    cols = ['chrom', 'start', 'end', 'sample']
    df.to_csv(out_file, columns=cols, sep="\t", header=False, index=False)

#-----------------------------------------------------------------------------------------------------------------------
# rm cal_accuracy fn here, & call instead from pipelines/rules/commons/evaluate_utils.py
#-----------------------------------------------------------------------------------------------------------------------
#def cal_accuracy(true_tracts, inferred_tracts):
#    """
#    Description:
#        Helper function for calculating accuracy.
#
#    Arguments:
#        true_tracts str: Name of the BED file containing true introgresssed tracts.
#        inferred_tracts str: Name of the BED file containing inferred introgressed tracts.
#
#    Returns:
#        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
#        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
#    """
#    truth_tracts = pybedtools.BedTool(true_tracts).sort().merge()
#    inferred_tracts =  pybedtools.BedTool(inferred_tracts).sort().merge()
#
#    total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
#    total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
#    true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])
#
#    if float(total_inferred_tracts) == 0: precision = np.nan
#    else: precision = true_positives / float(total_inferred_tracts) * 100
#    if float(total_true_tracts) == 0: recall = np.nan
#    else: recall = true_positives / float(total_true_tracts) * 100
#
#    return precision, recall

