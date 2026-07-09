# process pi per window
#cd /projects/psg/people/xcj768/arctic_wolves
#conda activate scikit-env
#---------------------------------------------------------------------------------------------------------------------
import pickle
import numpy as np
# read in pi objects per chrom, keep windows > 50,000 bases, output pi for these windows
	# [0]: pi per window, [1]: windows, [2]: n_bases per window, [3]: counts, number of variants per window
def proc_per_chrom(pref, chrom):
        inp_obj = open(f"genetic_diversity/nucl_diversity/{pref}_chr{chrom}_pi",'rb')
        chr_obj = pickle.load(inp_obj)
        # indices of windows to keep
        keep_index = chr_obj[2] > 50000
        keep_pi = chr_obj[0][keep_index]
        return keep_pi
#---------------------------------------------------------------------------------------------------------------------
# apply across all chromosomes 1-38, wolf only
pi_per_chrom = [proc_per_chrom('wolf_n56', chrom) for chrom in list(range(1,39))]
print([a.shape for a in pi_per_chrom])
#[(1227,), (854,), (919,), (883,), (889,), (776,), (810,), (743,), (611,), (693,), (744,), (725,), (632,), (610,), (642,), (596,), (643,), (558,), (537,), (581,), (509,), (614,), (523,), (477,), (516,), (390,), (459,), (412,), (418,), (402,), (399,), (388,), (314,), (421,), (265,), (308,), (309,), (239,)]
# flatten list of arrays, then take the mean of pi across windows
all_pi = np.concatenate(pi_per_chrom)
np.mean(all_pi)
#np.float64(0.0012858518286340856)
#---------------------------------------------------------------------------------------------------------------------
# full dataset
pref='n82.nooutgroup'
full_pi_per_chrom = [proc_per_chrom('n82.nooutgroup', chrom) for chrom in list(range(1,39))]
full_all_pi = np.concatenate(full_pi_per_chrom)
np.mean(full_all_pi)
#np.float64(0.0014403545382719873)
#---------------------------------------------------------------------------------------------------------------------
# estimate ne for singer
#π=4⋅Ne⋅m
#π/4m = Ne
#---------------------------------------------------------------------------------------------------------------------
mu = 4.96e-9
# wolves only dataset
np.mean(all_pi)/(4*mu) 
#np.float64(64811.08007228254)
# full dataset
np.mean(full_all_pi)/(4*mu)
#np.float64(72598.51503387032)
#---------------------------------------------------------------------------------------------------------------------
