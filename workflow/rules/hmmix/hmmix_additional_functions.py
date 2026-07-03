import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pybedtools


from helper_functions import sortby, Make_folder_if_not_exists, load_fasta, convert_to_bases, clean_files

from hmm_functions import HMMParam, get_default_HMM_parameters, write_HMM_to_file
from make_test_data import create_test_data
from hmm_functions import TrainModel, DecodeModel, HMMParam, read_HMM_parameters_from_file
from helper_functions import *
#from helper_functions import Load_observations_weights_mutrates


def Write_Decoded_output(outputprefix, segments, obs_file = None, admixpop_file = None, extrainfo = False):

    # Load archaic data
    if admixpop_file is not None:
        admix_pop_variants, admixpop_names = Annotate_with_ref_genome(admixpop_file, obs_file)

    # Are we doing haploid/diploid?
    outfile_mapper = {}
    for _, _, _, _, _, _, _, ploidity, _ in segments:
        if outputprefix == '/dev/stdout':
            outfile_mapper[ploidity] = '/dev/stdout'
        else:
            outfile_mapper[ploidity] = f'{outputprefix}.{ploidity}.txt'


    # Make output files and write headers
    outputfiles_handlers = defaultdict(str)
    for ploidity, output in outfile_mapper.items():

        Make_folder_if_not_exists(output)
        outputfiles_handlers[ploidity] = open(output, 'w')
        out = outputfiles_handlers[ploidity]

        if admixpop_file is not None:
            if extrainfo:
                out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tadmixpopvariants\t{}\tvariants\n'.format('\t'.join(admixpop_names)))
            else:
                out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tadmixpopvariants\t{}\n'.format('\t'.join(admixpop_names)))
        else:
            if extrainfo:
                out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\tvariants\n')
            else:
                out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\n')

    # Go through segments and write to output
    for chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, ploidity, variants in segments:

        out = outputfiles_handlers[ploidity]

        if admixpop_file is not None:
            archiac_variants_dict = defaultdict(int)
            for snp_position in variants.split(','):
                variant = admix_pop_variants[f'{chrom}_{snp_position}']
                if variant != '':
                    if '|' in variant:
                        for ind in variant.split('|'):
                            archiac_variants_dict[ind] += 1
                    else:
                        archiac_variants_dict[variant] += 1

                    archiac_variants_dict['total'] += 1

            archaic_variants = '\t'.join([str(archiac_variants_dict[x]) for x in ['total'] + admixpop_names])

            if extrainfo:
                print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, archaic_variants, variants, sep = '\t', file = out)
            else:
                print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, archaic_variants, sep = '\t', file = out)

        else:

            if extrainfo:
                print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, variants, sep = '\t', file = out)
            else:
                print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, sep = '\t', file = out)


    # Close output files
    for ploidity, out in outputfiles_handlers.items():
        out.close()


def sortby_haplotype(x):
    '''
    This function will sort haplotypes by number
    '''

    if '_hap' in x:
        return int(x.replace('_hap', ''))
    else:
        return x



def Load_observations_weights_mutrates(obs_file, weights_file, mutrates_file, window_size = 1000, haploid = False):

    obs_counter = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    haplotypes = defaultdict(int)

    with open(obs_file) as data:
        for line in data:
            if not line.startswith('chrom'):
                chrom, pos, ancestral_base, genotype = line.strip().split()
                rounded_pos = int(pos) - int(pos) % window_size

                if haploid:
                    for i, base in enumerate(genotype):
                        if base != ancestral_base:
                            obs_counter[chrom][rounded_pos][f'_hap{i+1}'].append(pos)
                            haplotypes[f'_hap{i+1}'] += 1
                else:
                    obs_counter[chrom][rounded_pos][''].append(pos)
                    haplotypes[''] += 1


    chroms, starts, variants, obs = [], [], [], []
    # In the case that there are NO derived variants - use weights to make list of zeros


    if len(obs_counter) == 0:

        if weights_file is None:
            #sys.exit(f'{obs_file} is empty! You need to provide a bed file!')
            raise Exception(f'{obs_file} is empty! You need to provide a bed file!')

        haplotypes[''] += 1
        callability = make_callability_from_bed(weights_file, window_size)
        for chrom in sorted(callability, key=sortby):
            lastwindow = max(callability[chrom]) + window_size

            for window in range(0, lastwindow, window_size):
                obs_counter[chrom][window][''].append('')
                chroms.append(f'{chrom}')
                starts.append(window)
                variants.append('')
                obs.append(0)


    # Otherwise fill out as normal
    else:
        for haplotype in sorted(haplotypes, key=sortby_haplotype):

            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    chroms.append(f'{chrom}{haplotype}')
                    starts.append(window)
                    variants.append(','.join(obs_counter[chrom][window][haplotype]))
                    obs.append(len(obs_counter[chrom][window][haplotype]))


    # Read weights file is it exists - else set all weights to 1
    if weights_file is None:
        weights = np.ones(len(obs))
    else:
        callability = make_callability_from_bed(weights_file, window_size)
        weights = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):

            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    weights.append(callability[chrom][window] / float(window_size))


    # Read mutation rate file is it exists - else set all mutation rates to 1
    if mutrates_file is None:
        mutrates = np.ones(len(obs))
    else:
        callability = make_callability_from_bed(mutrates_file, window_size)
        mutrates = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):

            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    mutrates.append(callability[chrom][window] / float(window_size))

    # Make sure there are no places with obs > 0 and 0 in mutation rate or weight
    for index, (observation, w, m) in enumerate(zip(obs, weights, mutrates)):
        if w*m == 0 and observation != 0:
            print(f'warning, you had {observation} observations but no called bases/no mutation rate at index:{index}. weights:{w}, mutrates:{m}')
            obs[index] = 0



    return np.array(obs).astype(int), chroms, starts, variants, np.array(mutrates).astype(float), np.array(weights).astype(float)



def train_hmm_individuals(tgt_ind_list, param, mutrates,  out_folder="", window_size=1000, haploid=False, weights=None):
    outfile_names = []
    infile_names = []
    for tgt_ind in tgt_ind_list:
        filename_start = "output_tgt_vs." + str(tgt_ind)
        infile = os.path.join(out_folder,  filename_start + ".txt")
        outfile = os.path.join(out_folder,  filename_start + "trained.json")
        train_hmm(infile, param, mutrates, outfile, window_size=1000, haploid=False, weights=None)
        infile_names.append(infile)
        outfile_names.append(outfile)
    return infile_names, outfile_names

def decode_hmm_individuals(comma_separated_tgt, trained_files,  mutrates, out_folder="", window_size=1000, haploid=False, weights=None):
    all_segments = []
    for outp_train_pair in zip(comma_separated_tgt, trained_files):
        print(outp_train_pair)
        filename_start = "output_tgt_vs." + str(outp_train_pair[0])
        infile = os.path.join(out_folder , filename_start +  ".txt")

        outfile = os.path.join(out_folder , filename_start + ".decoded.txt")
        #filename_start = "output_tgt_vs.txt." + str(tgt_ind)
        #train_hmm(filename_start + ".txt", param, mutrates, filename_start + "trained.json", window_size=1000, haploid=False, weights=None)
        new_df = decode_hmm(infile, outp_train_pair[1], mutrates, outfile, window_size=1000, haploid=False, weights=None)
        new_df["sample"] = outp_train_pair[0]
        all_segments.append(new_df)
    all_segments = pd.concat(all_segments)

    return all_segments


def train_hmm(obs, param, mutrates, out, window_size=1000, haploid=False, weights=None):
    try:
        hmm_parameters = read_HMM_parameters_from_file(param)
        obs, _, _, _, mutrates, weights = Load_observations_weights_mutrates(obs, weights, mutrates, window_size, haploid)

        print('-' * 40)
        print(hmm_parameters)
        print('> number of windows:', len(obs), '. Number of snps = ', sum(obs))
        print('> total callability:', round(np.sum(weights) / len(obs),2) )
        print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
        print('> Output is',out)
        print('> Window size is',window_size, 'bp')
        print('> Haploid',haploid)
        print('-' * 40)

        hmm_parameters = TrainModel(obs, mutrates, weights, hmm_parameters)
        write_HMM_to_file(hmm_parameters, out)
    except:
        print("training of model failed - probably due to no derived observations!")



def decode_hmm(obs, param, mutrates, out, window_size=1000, haploid=False, weights=None):
        try:
            obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(obs, weights, mutrates, window_size, haploid)
            hmm_parameters = read_HMM_parameters_from_file(param)

            print('-' * 40)
            print(hmm_parameters)
            print('> number of windows:', len(obs), '. Number of snps = ', sum(obs))
            print('> total callability:', round(np.sum(weights) / len(obs),2) )
            print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
            print('> Output prefix is',out)
            print('> Window size is',window_size, 'bp')
            print('> Haploid',haploid)
            print('-' * 40)

            # Find segments and write output


            segments = DecodeModel(obs, chroms, starts, variants, mutrates, weights, hmm_parameters)
            #Write_Decoded_output(out, segments, obs, admixpop, extrainfo)

            Write_Decoded_output(out, segments, obs, None, False)
            if haploid:
                segments_df = pd.read_table(out + ".haploid.txt")
            else:
                segments_df = pd.read_table(out + ".diploid.txt")
            print(segments_df)


            #return probabilities
            return segments_df
        except:
            print("decoding of model failed - probably due to no derived observations!")
            segments_df = pd.DataFrame()
            return segments_df


def convert_to_bases_set(genotype, both_bases, ref_set= ["0", "1"]):

    return_genotype = 'NN'
    separator = None

    if '/' in genotype or '|' in genotype:
        separator = '|' if '|' in genotype else '/'

        base1, base2 = [x for x in genotype.split(separator)]
        if base1.isnumeric() and base2.isnumeric():
            base1, base2 = int(base1), int(base2)

            if both_bases[base1] in ref_set and both_bases[base2] in ref_set:
                return_genotype = both_bases[base1] + both_bases[base2]

    return return_genotype


def make_ingroup_custom_mut(ingroup_individuals, bedfile, vcffiles, outprefix, outgroupfile, ancestralfiles, ref_set=["0", "1"]):

    # handle output files
    Make_folder_if_not_exists(outprefix)
    outfile_handler = defaultdict(str)
    for individual in ingroup_individuals:
        outfile_handler[individual] = open(f'{outprefix}.{individual}.txt','w')
        print('chrom', 'pos', 'ancestral_base', 'genotype', sep = '\t', file = outfile_handler[individual])

    individuals_for_bcf = ','.join(ingroup_individuals)

    for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):

        if ancestralfile is not None:
            ancestral_allele = load_fasta(ancestralfile)

        if bedfile is not None:
            command = f'bcftools view -v snps -s {individuals_for_bcf} -T {bedfile} {vcffile} | bcftools norm -m +any | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'
        else:
            command = f'bcftools view -v snps -s {individuals_for_bcf} {vcffile} | bcftools norm -m +any | vcftools --vcf - --exclude-positions {outgroupfile} --recode --stdout'

        print('Running command:')
        print(command, '\n\n')

        for index, line in enumerate(os.popen(command)):

            if line.startswith('#CHROM'):
                individuals_in_vcffile = line.strip().split()[9:]

            #print("show the lines")
            #print(line)

            if not line.startswith('#'):

                chrom, pos, _, ref_allele, alt_allele = line.strip().split()[0:5]
                pos = int(pos)
                genotypes = [x.split(':')[0] for x in line.strip().split()[9:]]
                all_bases = [ref_allele] + alt_allele.split(',')

                if ref_allele in ref_set:

                    for original_genotype, individual in zip(genotypes, individuals_in_vcffile):
                        genotype = convert_to_bases_set(original_genotype, all_bases, ref_set)


                        if ancestralfile is not None:
                            # With ancestral information look for derived alleles
                            ancestral_base = ancestral_allele[pos-1]
                            if ancestral_base in all_bases and genotype.count(ancestral_base) != 2 and genotype != 'NN':
                                print(chrom, pos, ancestral_base, genotype, sep = '\t', file = outfile_handler[individual])

                        else:
                            # If no ancestral information is provided only include heterozygous variants
                            if genotype[0] != genotype[1]:
                                print(chrom, pos, ref_allele, genotype, sep = '\t', file = outfile_handler[individual])


                if index % 100000 == 0:
                    print(f'at line {index} at chrom {chrom} and position {pos}')

    # Clean log files generated by vcf and bcf tools
    clean_files('out.log')

    for individual in ingroup_individuals:
        outfile_handler[individual].close()


def make_out_group_custom_mut(individuals_input, bedfile, vcffiles, outputfile, ancestralfiles, refgenomefiles, ref_set=["0", "1"]):

    Make_folder_if_not_exists(outputfile)
    outgroup_individuals = ','.join(individuals_input)

    with open(outputfile + '.unsorted', 'w') as out:

        print('chrom', 'pos', 'ref_allele_info', 'alt_allele_info', 'ancestral_base', sep = '\t', file = out)

        for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):

            if ancestralfile is not None:
                ancestral_allele = load_fasta(ancestralfile)

            if bedfile is not None:
                command = f'bcftools view -s {outgroup_individuals} -T {bedfile} {vcffile} | bcftools norm -m -any | bcftools view -v snps | vcftools --vcf - --counts --stdout'
            else:
                command = f'bcftools view -s {outgroup_individuals} {vcffile} | bcftools norm -m -any | bcftools view -v snps | vcftools --vcf - --counts --stdout'

            print(f'Processing {vcffile}...')
            print('Running command:')
            print(command, '\n\n')

            variants_seen = defaultdict(int)
            for index, line in enumerate(os.popen(command)):
                if not line.startswith('CHROM'):

                    chrom, pos, _, _, ref_allele_info, alt_allele_info = line.strip().split()

                    ref_allele, ref_count = ref_allele_info.split(':')
                    alt_allele, alt_count = alt_allele_info.split(':')
                    pos, ref_count, alt_count = int(pos),  int(ref_count), int(alt_count)

                    # Always include polymorphic sites
                    if alt_count * ref_count > 0:
                        ancestral_base = ref_allele if ref_count > alt_count else alt_allele

                        # Use ancestral base info if available
                        if ancestralfile is not None:
                            ancestral_base_temp = ancestral_allele[pos-1]
                            if ancestral_base_temp in [ref_allele, alt_allele]:
                                 ancestral_base = ancestral_base_temp

                        print(chrom, pos, ref_allele_info, alt_allele_info, ancestral_base, sep = '\t', file = out)
                        variants_seen[pos-1] = 1

                    # Fixed sites
                    elif alt_count * ref_count == 0:
                        ancestral_base = ref_allele if ref_count > alt_count else alt_allele

                        # Use ancestral base info if available
                        if ancestralfile is not None:
                            ancestral_base_temp = ancestral_allele[pos-1]
                            if ancestral_base_temp in [ref_allele, alt_allele]:
                                 ancestral_base = ancestral_base_temp

                        if ancestral_base == alt_allele:
                            derived_count = ref_count
                        else:
                             derived_count = alt_count

                        if derived_count > 0:
                            print(chrom, pos, ref_allele_info, alt_allele_info, ancestral_base, sep = '\t', file = out)
                            variants_seen[pos-1] = 1


                    if index % 100000 == 0:
                        print(f'at line {index} at chrom {chrom} and position {pos}')

            # If reference genome is provided then remove positions where the reference and ancestral differ AND which is not found in the outgroup
            if reffile is not None and ancestralfile is not None:
                print('Find fixed derived sites')
                refgenome_allele = load_fasta(reffile)

                for index, (refbase, ancbase) in enumerate(zip(refgenome_allele, ancestral_allele)):
                    if ancbase in ref_set  and refbase in ref_set :
                        if refbase != ancbase and variants_seen[index] == 0:
                            print(chrom, index + 1, f'{refbase}:100', f'{ancbase}:0', ancbase, sep = '\t', file = out)

    # Sort outgroup file
    print('Sorting outgroup file')
    positions_to_sort = defaultdict(lambda: defaultdict(str))
    with open(outputfile + '.unsorted') as data, open(outputfile, 'w') as out:
        for line in data:
            if line.startswith('chrom'):
                out.write(line)
            else:
                chrom, pos = line.strip().split()[0:2]
                positions_to_sort[chrom][int(pos)] = line

        for chrom in sorted(positions_to_sort, key=sortby):
            for pos in sorted(positions_to_sort[chrom]):
                line =  positions_to_sort[chrom][pos]
                out.write(line)

    # Clean log files generated by vcf and bcf tools
    clean_files(outputfile + '.unsorted')
    clean_files('out.log')




def cal_accuracy(true_tracts, inferred_tracts):
    """
    Description:
        Helper function for calculating accuracy.

    Arguments:
        true_tracts str: Name of the BED file containing true introgresssed tracts.
        inferred_tracts str: Name of the BED file containing inferred introgressed tracts.

    Returns:
        precision float: Amount of true introgressed tracts detected divided by amount of inferred introgressed tracts.
        recall float: Amount ot true introgressed tracts detected divided by amount of true introgressed tracts.
    """
    print("show files")
    print(true_tracts)
    print(inferred_tracts)


    try:
        truth_tracts = pybedtools.BedTool(true_tracts).sort().merge()
        inferred_tracts =  pybedtools.BedTool(inferred_tracts).sort().merge()

        total_inferred_tracts = sum([x.stop - x.start for x in (inferred_tracts)])
        total_true_tracts =  sum([x.stop - x.start for x in (truth_tracts)])
        true_positives = sum([x.stop - x.start for x in inferred_tracts.intersect(truth_tracts)])

        if float(total_inferred_tracts) == 0: precision = np.nan
        else: precision = true_positives / float(total_inferred_tracts) * 100
        if float(total_true_tracts) == 0: recall = np.nan
        else: recall = true_positives / float(total_true_tracts) * 100

        return precision, recall

    except:
        return 0, 0



def process_output(segments_df, outfolder, out_file_prefix, src_id, cutoff_list=[0.01, 0.25, 0.5, 0.75, 0.99], return_filenames=False):
    segments_df.loc[segments_df['state'] != src_id, 'mean_prob'] = 1 - segments_df.loc[segments_df['state'] != src_id, 'mean_prob']

    if return_filenames:
        filenames = []

    for cutoff in cutoff_list:
        cutoff_df = segments_df.copy()
        cutoff_df = cutoff_df.drop('state', axis=1)
        cutoff_df = cutoff_df[cutoff_df['mean_prob'] > cutoff]

        cols = ['chrom', 'start', 'end', 'sample']

        new_file = os.path.join(outfolder, out_file_prefix + str(cutoff) + ".bed")
        cutoff_df.to_csv(new_file, columns=cols, sep="\t", header=False, index=False)

        if return_filenames:
            filenames.append([cutoff, new_file])

    if return_filenames:
        return filenames
