import argparse
from scipy.stats import norm
from scipy.stats import nbinom

parser = argparse.ArgumentParser(description="input params from smk")
parser.add_argument("--seqlen", type=float, help="seq_len")
parser.add_argument("--mutrate", type=float, help="mut_rate")
parser.add_argument("--recrate", type=float, help="rec_rate")
parser.add_argument("--outfile", type=str, help="outfile")

args = parser.parse_args()

seq_len = args.seqlen
N0 = 1000
mut_rate_mean = args.mutrate
rec_rate_mean = args.recrate
output = args.outfile

scaled_mut_rate_mean = 4*N0*mut_rate_mean*seq_len
scaled_mut_rate_sdv = 0.233

scaled_rec_rate_mean = 4*N0*rec_rate_mean*seq_len
mut_rate_list = norm.rvs(loc=scaled_mut_rate_mean, scale=scaled_mut_rate_sdv, size=20000)
rec_rate_list = nbinom.rvs(n=0.5, p=0.5/(0.5+scaled_rec_rate_mean), size=20000)

with open(output, 'w') as o:
    for i in range(len(mut_rate_list)):
        if mut_rate_list[i] < 0.001: mut_rate_list[i] = 0.001
        if rec_rate_list[i] < 0.001: rec_rate_list[i] = 0.001
        mut_rate_new = mut_rate_list[i]
        rec_rate_new = rec_rate_list[i]
        o.write(f'{mut_rate_new}\t{rec_rate_new}\n')