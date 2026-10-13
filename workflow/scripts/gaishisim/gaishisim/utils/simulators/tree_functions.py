import allel
import tskit
import numpy as np


def compute_diversity(vcf_file, chr="1"):
    callset = allel.read_vcf(vcf_file, region=chr, fields=['CHROM','POS','GT'])
    pos = callset['variants/POS']
    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()
    # calc nucleotide diversity in windows, size 100kb

    pis = []
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=100000)
    pis.append(pi)
    final_pi = np.mean(pis)
    return final_pi




def compute_eff_popsize_from_diversity(mutation_rate, vcf_file=None, ts=None, chr="1"):


    if ts:
        pi = ts.diversity()
    elif vcf_file:
        pi = compute_diversity(vcf_file, chr="1")
    else:
        raise Exception()

    neff = pi / (4 * mutation_rate)
    return neff

