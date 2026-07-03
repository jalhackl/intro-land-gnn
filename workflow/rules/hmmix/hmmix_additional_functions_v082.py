import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pybedtools
import sys
import pickle

# path to skov helper scripts for hmmix
skov_dir = "resources/skov_scripts/Introgression-detection/src/"
sys.path.append(skov_dir)


#fn DecodeModel -> replaced w Calculate_Posterior_probabillities
from hmm_functions import TrainModel, write_HMM_to_file, read_HMM_parameters_from_file, Write_Decoded_output, Calculate_Posterior_probabillities, PMAP_path, Viterbi_path, Hybrid_path, Convert_genome_coordinates, Write_posterior_probs, Make_inhomogeneous_transition_matrix, Simulate_from_transition_matrix, Write_inhomogeneous_transition_matrix, Emission_probs_poisson

# create_test_data -> simulate_path, write_data
from make_test_data import simulate_path, write_data
from helper_functions import *

# check decode_hmm function *

def Write_posterior_probs(chroms, starts, weights, mutrates, post_seq, path, variants, hmm_parameters, filename):
    post_seq = post_seq.T
    with open(filename, 'w') as out:
        state_names = '\t'.join(hmm_parameters.state_names)
        print('chrom', 'start', 'called_sequence', 'mutationrate', state_names,'state','variants', sep = '\t', file = out)
        for (chrom, start, w, m, posterior, state, var) in zip(chroms, starts, weights, mutrates, post_seq, path, variants):
            posterior_to_print = '\t'.join([str(round(x, 4)) for x in posterior])
            print(chrom, start, w, m, posterior_to_print, hmm_parameters.state_names[state], var, sep = '\t', file = out)


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
            raise ValueError(f'{obs_file} is empty! You need to provide a bed file!')
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
        filename_start = f"output_tgt_vs.{str(tgt_ind)}"
        infile = os.path.join(out_folder,  filename_start + ".txt")
        outfile = os.path.join(out_folder,  filename_start + "trained.json")
        train_hmm(infile, param, mutrates, outfile, window_size=1000, haploid=False, weights=None)
        infile_names.append(infile)
        outfile_names.append(outfile)
    return infile_names, outfile_names


#not used, but could be advantageous in case one wants to create all ind-files even in case of failing
def train_hmm_individuals_with_fix(tgt_ind_list, param, mutrates,
                          out_folder="", window_size=1000,
                          haploid=False, weights=None):

    outfile_names = []
    infile_names = []

    for tgt_ind in tgt_ind_list:

        filename_start = f"output_tgt_vs.{str(tgt_ind)}"
        infile = os.path.join(out_folder, filename_start + ".txt")
        outfile = os.path.join(out_folder, filename_start + "trained.json")
        skipfile = os.path.join(out_folder, filename_start + ".skip")

        # --------------------------------------------------
        # CASE 1: missing or empty observation file
        # --------------------------------------------------
        if (not os.path.exists(infile)) or os.path.getsize(infile) == 0:

            print(f"[SKIP] {tgt_ind}: no SNPs → writing skip model")

            import json
            with open(outfile, "w") as f:
                json.dump({"status": "skip", "individual": tgt_ind}, f)

            open(skipfile, "w").write("skip\n")

            infile_names.append(infile)
            outfile_names.append(outfile)
            continue

        # --------------------------------------------------
        # CASE 2: normal training
        # --------------------------------------------------
        try:
            train_hmm(infile, param, mutrates, outfile,
                      window_size=window_size,
                      haploid=haploid,
                      weights=weights)

        except Exception as e:
            print(f"[WARN] training failed for {tgt_ind}: {e}")

            import json
            with open(outfile, "w") as f:
                json.dump({"status": "skip", "reason": "training_failed"}, f)

            open(skipfile, "w").write("skip\n")

        infile_names.append(infile)
        outfile_names.append(outfile)

    return infile_names, outfile_names


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

        hmm_parameters = TrainModel(obs, mutrates, weights, hmm_parameters, epsilon=1e-3, maxiterations=1000)
        write_HMM_to_file(hmm_parameters, out)
    except Exception as e:
        print("training of model failed -", e)



#-----------------------------------------------------------------------------------------------------------------------
# if  hybrid=-1, viterbi=False => use default posterior decoding algorithm
    # contains pickle step to output intermediate output
    # should take one or multiple chroms
def decode_hmm(obs, param, mutrates, out, chrom, window_size=1000, haploid=False, weights=None, hybrid=-1, viterbi=False):
        try:
            obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(obs, weights, mutrates, window_size, haploid)
            hmm_parameters = read_HMM_parameters_from_file(param)
            CHROMOSOME_BREAKPOINTS = [x for x in find_runs(chroms)]
            print('-' * 40)
            print(hmm_parameters)
            print(f'> chromosomes to use: {chrom}')
            print('> number of windows:', len(obs), '. Number of snps = ', sum(obs))
            print('> total callability:', round(np.sum(weights) / len(obs),2) )
            print('> average mutation rate per bin:', round(np.sum(mutrates * weights) / np.sum(weights), 2) )
            print('> Output prefix is',out)
            print('> Window size is',window_size, 'bp')
            print('> Haploid',haploid)
            print('-' * 40)
            # find segments
            # emissions, posterior prob
            emissions = Emission_probs_poisson(hmm_parameters.emissions, obs, weights, mutrates)
            posterior_probs = Calculate_Posterior_probabillities(emissions, hmm_parameters)
            print('check emissions',emissions)
            print('check posteriorprob',posterior_probs)
            # decoding
            if hybrid != -1:
                if 0 <= hybrid <= 1:
                    print(f'> Decode using hybrid algorithm with parameter: {hybrid}')
                    print('-' * 40)
                    logged_posterior_probs = np.log(posterior_probs.T)
                    path = Hybrid_path(emissions, hmm_parameters.starting_probabilities, hmm_parameters.transitions, logged_posterior_probs, hybrid)
                else:
                    sys.exit('\n\nERROR! Hybrid parameter must be between 0 and 1\n\n')
            else:
                if viterbi:
                    print('> Decode using viterbi algorithm')
                    print('-' * 40)
                    path = Viterbi_path(emissions, hmm_parameters)
                else:
                    print('> Decode with posterior decoding (PMAP)')
                    print('-' * 40)
                    path = PMAP_path(posterior_probs)
                    # save posterior prob
            Write_posterior_probs(chroms, starts, weights, mutrates, posterior_probs, path, variants, hmm_parameters, out)
            segments = Convert_genome_coordinates(window_size, CHROMOSOME_BREAKPOINTS, starts, variants, posterior_probs, path, hmm_parameters, weights, mutrates, obs)
            print('check path',path)
            print('check segments', segments)
            with open(f"{out}.test", "wb") as f:
                pickle.dump((emissions, posterior_probs, segments, obs), f)
            test_Write_Decoded_output(out, segments, obs)
            if haploid:
                segments_df = pd.read_table(out + ".haploid.txt")
            else:
                segments_df = pd.read_table(out + ".diploid.txt")
            print('check segments df',segments_df)
            return segments_df
        except Exception as e:
            print("decoding of model failed -", e)


#-----------------------------------------------------------------------------------------------------------------------

def decode_hmm_individuals(comma_separated_tgt, trained_files, mutrates, chrom,
                           out_folder="", window_size=1000,
                           haploid=False, weights=None,
                           hybrid=-1, viterbi=False):

    all_segments = []

    for ind, trained_file in zip(comma_separated_tgt, trained_files):

        print(f"Decoding individual: {ind}")

        filename_start = f"output_tgt_vs.{str(ind)}"
        infile = os.path.join(out_folder, filename_start + ".txt")
        outfile = os.path.join(out_folder, filename_start + ".decoded.txt")

        new_df = decode_hmm(
            infile,
            trained_file,
            mutrates,
            outfile,
            chrom,
            window_size=window_size,
            haploid=haploid,
            weights=weights,
            hybrid=hybrid,
            viterbi=viterbi
        )

        # Skip individuals with failed decoding (e.g. zero SNPs)
        if new_df is None:
            print(f"Skipping {ind} — no informative SNPs found in .txt-file.")
            continue

        new_df["sample"] = ind
        all_segments.append(new_df)

    # If all individuals failed, return empty DataFrame
    if len(all_segments) == 0:
        print("WARNING: No individuals could be decoded.")
        return pd.DataFrame()

    return pd.concat(all_segments, ignore_index=True)

#-----------------------------------------------------------------------------------------------------------------------
def artemis_alpha(param, out, out_plot, windows=500000, iterations=10, start=0, end=1, steps=101, seed=42):
            hmm_parameters = read_HMM_parameters_from_file(param)
            print('-' * 40)
            print(hmm_parameters)  
            print(f'> Save data to {out}')
            print(f'> Save plot to {out_plot}')
            print(f'> Number of windows: {windows}')
            print(f'> Number of iterations: {iterations}')
            print(f'> Test {steps} alphas between {start} and {end}:')
            print(f'> Seed is {seed}')
            print('-' * 40)
            out_alpha = Find_best_alpha(hmm_parameters, windows, out, out_plot, iterations,  start, end, steps, seed)
            return out_alpha
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def test_Write_Decoded_output(outputprefix, segments, obs_file = None):
    try:
        # check ploidy
        outfile_mapper = {}
        for _, _, _, _, _, _, _, ploidity, _, _, _ in segments:
            if outputprefix == '/dev/stdout':
                outfile_mapper[ploidity] = '/dev/stdout'
            else:
                outfile_mapper[ploidity] = f'{outputprefix}.{ploidity}.txt'
        # Create output files
        outputfiles_handlers = {}
        for ploidity, output in outfile_mapper.items():
            Make_folder_if_not_exists(output)
            outputfiles_handlers[ploidity] = open(output, 'w')
            out = outputfiles_handlers[ploidity]
            # Write headers
            out.write('chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\n')
        # Process and write each segment to the appropriate output file
        for chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, ploidity, _, _, _ in segments:
            out = outputfiles_handlers[ploidity]
            print(chrom, genome_start, genome_end, genome_length, state, mean_prob, snp_counter, sep='\t', file=out)
    except Exception as e:
        print(f"writing decoded output failed: {e}")    
    finally:
        # Close all output file handlers
        for out in outputfiles_handlers.values():
            out.close()

def Make_folder_if_not_exists(path):
    '''
    Check if path exists - otherwise make it;
    '''
    path = os.path.dirname(path)
    if path != '':
        if not os.path.exists(path):
            os.makedirs(path)


#-----------------------------------------------------------------------------------------------------------------------


def convert_to_bases_set(genotype, both_bases, ref_set=None):
    if ref_set is None:
        ref_set = ["0", "1"]
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

#-----------------------------------------------------------------------------------------------------------------------


def process_output(segments_df, outfolder, out_file_prefix, src_id, cutoff_list=None, return_filenames=False):
    segments_df.loc[segments_df['state'] != src_id, 'mean_prob'] = 1 - segments_df.loc[segments_df['state'] != src_id, 'mean_prob']
    if cutoff_list is None:
        cutoff_list = [0.01, 0.25, 0.5, 0.75, 0.99]
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
