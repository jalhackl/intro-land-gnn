#conda activate scikit-env
#cd /projects/psg/people/xcj768/arctic_wolves
# 1) calculate proportion of callable sites per 100kb window [windows used to calculate nucleotide diversity]
#---------------------------------------------------------------------------------------------------------------------
import pyranges as pr
import pickle
import numpy as np
import pandas as pd
#---------------------------------------------------------------------------------------------------------------------
# good sites from the called gts dataset
goodsites_path='site_filters/good_sites_depth_hete0_norep.bed'
goodsites_gr = pr.read_bed(goodsites_path)
#goodsites_gr
#index    |    Chromosome    Start     End       Name      Score    Strand
#int64    |    category      int64     int64     object    int64    category
#-------  ---  ------------  --------  --------  --------  -------  ----------
#0        |    chr1          212105    217661    .         0        .
#1        |    chr1          217664    217667    .         0        .
#2        |    chr1          217734    218015    .         0        .
#3        |    chr1          218016    218021    .         0        .
#...      |    ...           ...       ...       ...       ...      ...
#2301508  |    chr38         23912969  23914398  .         0        .
#2301509  |    chr38         23914399  23914401  .         0        .
#2301510  |    chr38         23914402  23914404  .         0        .
#2301511  |    chr38         23914405  23914408  .         0        .
#PyRanges with 2301512 rows, 6 columns, and 1 index columns.
#Contains 38 chromosomes and 1 strands (including non-genomic strands: .).
#---------------------------------------------------------------------------------------------------------------------
# read in pi objects per chrom
       # [0]: pi per window, [1]: windows, [2]: n_bases per window, [3]: counts, number of variants per window
def proc_per_chrom(pref, chrom):
        inp_obj = open(f"genetic_diversity/nucl_diversity/{pref}_chr{chrom}_pi",'rb')
        chr_obj = pickle.load(inp_obj)
        return chr_obj
#---------------------------------------------------------------------------------------------------------------------
# apply across all chromosomes 1-38, wolf only
#pi_per_chrom = [proc_per_chrom('wolf_n56', chrom) for chrom in list(range(1,39))]
# or for full dataset -
pi_per_chrom = [proc_per_chrom('n82.nooutgroup', chrom) for chrom in list(range(1,39))]
#---------------------------------------------------------------------------------------------------------------------
# testing for chr1, ie 0 = chrom1, 1 = windows
#>>> pi_per_chrom[0][1]
#array([[     1763,    101762],
#       [   101763,    201762],
#       [   201763,    301762],
#       ...,
#       [122401763, 122501762],
#       [122501763, 122601762],
#       [122601763, 122678715]], shape=(1227, 2))
#---------------------------------------------------------------------------------------------------------------------
# convert windows for one chr to pandas df then to pyranges
# overlap with the good sites  
# calc proportion of windows overlapping good sites
# output informative windows per chrom
#---------------------------------------------------------------------------------------------------------------------
# generate for all windows for this chrom
#---------------------------------------------------------------------------------------------------------------------
def inform_fraction_per_window(int_wind, orig_winds, uniquewinds, inp):
        # extract all intersect ranges for 1 original window
        s1 = int_wind[int_wind['orig_win'] == uniquewinds.iloc[inp]]
        # extract coords of original window
        orig_coords = orig_winds[orig_winds['win_index'] == uniquewinds.iloc[inp]]
        #proportion of window covered: sum of intersect regions / length of window
        orig_coords['proportion'] = sum(s1['lengths'])/100000
        return orig_coords
#---------------------------------------------------------------------------------------------------------------------
def proc_one_chrom(pref, inpchrom):
        # extract single chrom, windows for which pi calculated + covert to df
        chrom1_df = pd.DataFrame({'Start': pi_per_chrom[inpchrom][1][:, 0], 'End': pi_per_chrom[inpchrom][1][:, 1]})
        # add chr identifier
        chrom1_df.insert(0, 'Chromosome', f"chr{inpchrom+1}")
        # convert to pyranges
        chr1_wind=pr.PyRanges(chrom1_df)
        chr1_wind['win_index'] = chr1_wind.index
        # intersect chr1 w good sites for chr1
        ints = chr1_wind.intersect_overlaps(goodsites_gr)
        # add annotations - widths per window, which window pertains to & window coords
        ints['lengths'] = ints.lengths()
        ints['orig_win'] = ints.index
        # extract unique window identifiers
        unique_winds = ints["orig_win"].drop_duplicates()
        # run fn inform_fraction_per_window for all windows for this chromosome
        unique_windows_chrom1 = [inform_fraction_per_window(ints, chr1_wind, unique_winds, i) for i in list(range(0,len(unique_winds)))]
        all_frac_df = pd.concat(unique_windows_chrom1)
        # original window indices that did not have any overlap with the callable sites
        no_overlap_coords = chr1_wind[~chr1_wind['win_index'].isin(unique_winds)]
        # set proportion to 0
        no_overlap_coords = no_overlap_coords.assign(proportion=0)
        # concatenate these 2 dfs
        both_dfs = pd.concat([all_frac_df, no_overlap_coords])
        # sort by window index
        both_sorted = both_dfs.sort_values(by='win_index')
        with open(f"genetic_diversity/nucl_diversity/callable_info_per_window/{pref}_{inpchrom+1}", "wb") as f:
              pickle.dump((both_sorted), f)
        return both_sorted
#---------------------------------------------------------------------------------------------------------------------
#all_chrom_windows = [proc_one_chrom('wolf_n56', chrom) for chrom in list(range(0,38))]
# generate also for 'n82.nooutgroup'
all_chrom_windows = [proc_one_chrom('n82.nooutgroup', chrom) for chrom in list(range(0,38))]
#---------------------------------------------------------------------------------------------------------------------
# same windows for both datasets
#>>> all_chrom_windows[0]
#index    |    Chromosome    Start      End        win_index    proportion
#int64    |    object        int64      int64      int64        float64
#-------  ---  ------------  ---------  ---------  -----------  ------------
#0        |    chr1          1763       101762     0            0.0
#1        |    chr1          101763     201762     1            0.0
#2        |    chr1          201763     301762     2            0.85626
#3        |    chr1          301763     401762     3            0.94636
#...      |    ...           ...        ...        ...          ...
#1223     |    chr1          122301763  122401762  1223         0.87688
#1224     |    chr1          122401763  122501762  1224         0.86359
#1225     |    chr1          122501763  122601762  1225         0.94771
#1226     |    chr1          122601763  122678715  1226         0.71469
#PyRanges with 1227 rows, 5 columns, and 1 index columns.
#Contains 1 chromosomes.
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# estimate ne for singer
#π=4⋅Ne⋅m
#π/4m = Ne
mu = 4.96e-9
#---------------------------------------------------------------------------------------------------------------------
# per chrom filter pi by windows with sufficient proportion of callable sites
def filter_pi_by_windows(inpchrom):
        wins_keep = all_chrom_windows[inpchrom].loc[all_chrom_windows[inpchrom]['proportion'] > 0.7, 'win_index']
        chrom_pi =  pi_per_chrom[inpchrom][0]
        # subset pi by windows with sufficient proportion of callable sites
        keep_pi = chrom_pi[wins_keep]
        return keep_pi
#---------------------------------------------------------------------------------------------------------------------
# full dataset, filter by proportion > 0.7
full_pi_keep = list(map(filter_pi_by_windows, list(range(0,38))))
full_all_pi = np.concatenate(full_pi_keep)
np.mean(full_all_pi)
#np.float64(0.0014352918319049694) # pi
np.mean(full_all_pi)/(4*mu)
#np.float64(72343.33830166176) # ne
#---------------------------------------------------------------------------------------------------------------------
# try more stringent proportion for fraction of callable sites 
def filter_pi_by_windows(inpchrom):
        wins_keep = all_chrom_windows[inpchrom].loc[all_chrom_windows[inpchrom]['proportion'] > 0.9, 'win_index']
        chrom_pi =  pi_per_chrom[inpchrom][0]
        # subset pi by windows with sufficient proportion of callable sites
        keep_pi = chrom_pi[wins_keep]
        return keep_pi
#---------------------------------------------------------------------------------------------------------------------
# full dataset, filter by proportion > 0.9
full_pi_keep = list(map(filter_pi_by_windows, list(range(0,38))))
full_all_pi = np.concatenate(full_pi_keep)
np.mean(full_all_pi)
#np.float64(0.0013862145738430915) #pi
np.mean(full_all_pi)/(4*mu)
#np.float64(69869.68618160744) # ne
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# wolf only dataset, filter by proportion > 0.7
full_pi_keep = list(map(filter_pi_by_windows, list(range(0,38))))
full_all_pi = np.concatenate(full_pi_keep)
np.mean(full_all_pi)
#np.float64(0.0012796586851183665) # pi
np.mean(full_all_pi)/(4*mu)
#np.float64(64498.925661207984) # ne
#---------------------------------------------------------------------------------------------------------------------
# wolf only dataset, filter by proportion > 0.9
full_pi_keep = list(map(filter_pi_by_windows, list(range(0,38))))
full_all_pi = np.concatenate(full_pi_keep)
np.mean(full_all_pi)
# np.float64(0.0012364175206735848) # pi
np.mean(full_all_pi)/(4*mu)
#np.float64(62319.43148556375) # ne
#---------------------------------------------------------------------------------------------------------------------

