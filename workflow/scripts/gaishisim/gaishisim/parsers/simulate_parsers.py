# Copyright 2024 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import argparse
import sys
from gaishisim.parsers.argument_validation import positive_int
from gaishisim.parsers.argument_validation import positive_number
from gaishisim.parsers.argument_validation import existed_file


def _run_simulation(args: argparse.Namespace) -> None:
    """
    Executes a simulation process with specified parameters.

    Parameters
    ----------
    args : argparse.Namespace
        A namespace object obtained from argparse, containing simulation parameters:
        - demes: File path to the demographic model specification in YAML format.
        - replicate: Number of simulation replicates.
        - nref: Size of the reference population.
        - ntgt: Size of the target population.
        - ref_id: Identifier for the reference population.
        - tgt_id: Identifier for the target population.
        - src_id: Identifier for the source population.
        - ploidy: Ploidy of the organisms being simulated.
        - is_phased: Indicates if the simulated data should be phased.
        - seq_len: Length of the sequence to simulate.
        - mut_rate: Mutation rate to use in the simulation.
        - rec_rate: Recombination rate to use in the simulation.
        - nprocess: Number of processes to use for parallel simulations.
        - output_prefix: Prefix for output files.
        - output_dir: Directory where output files will be saved.
        - seed: Random seed for reproducibility.
        - nfeature: Number of features to simulate.
        - num_polymorphisms: Number of polymorphisms in each genotype matrix to simulate.
        - num_upsamples: Number of samples after upsampling.
        - output_h5: Boolean flag to save output in HDF5 format.
        - is_sorted: Boolean flag to indicate whether to sort the genotype matrices.
        - only_intro: Boolean flag to simulate only introgressed fragments.
        - only_non_intro: Boolean flag to simulate only non-introgressed fragments.
        - force_balanced: Boolean flag to ensure a balanced distribution of introgressed and
                          non-introgressed classes in the training data.
        - keep_sim_data: Boolean flag to keep or discard simulation data.
        - chunk_size: integer for HDF chunk size
        - mut_map
        - rec_map
        - start_rep

    """

    import demes
    from gaishisim.simulate import simulate_test_data

    demog = demes.load(args.demes)
    pops = [d.name for d in demog.demes]
    if args.ref_id not in pops:
        print(
            f"gaishisim simulate: error: argument --ref_id: Population {args.ref_id} is not found in the demographic model file {args.demes}"
        )
        sys.exit(1)
    if args.tgt_id not in pops:
        print(
            f"gaishisim simulate: error: argument --tgt_id: Population {args.tgt_id} is not found in the demographic model file {args.demes}"
        )
        sys.exit(1)
    if args.src_id not in pops:
        print(
            f"gaishisim simulate: error: argument --src_id: Population {args.src_id} is not found in the demographic model file {args.demes}"
        )
        sys.exit(1)

    simulate_test_data(
        demo_model_file=args.demes,
        nrep=args.replicate,
        nref=args.nref,
        ntgt=args.ntgt,
        ref_id=args.ref_id,
        tgt_id=args.tgt_id,
        src_id=args.src_id,
        ploidy=args.ploidy,
        is_phased=args.phased,
        seq_len=args.seq_len,
        mut_rate=args.mut_rate,
        rec_rate=args.rec_rate,
        nprocess=args.nprocess,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        seed=args.seed,
        mutation_model=args.mutation_model,
        rec_map=args.rec_map,
        mut_map=args.mut_map,
        L_mut_map = args.L_mut_map,
        cv_mut_map = args.cv_mut_map,
        shuffle_rec_map=args.shuffle_rec_map,
        true_tracts_batch_processing=args.true_tracts_batch_processing,
        # in case one wants to start with a specific subfolder nr.
        start_rep=args.start_rep
    )


def add_simulate_parsers(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for simultating data.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.

    """

    parser = subparsers.add_parser(
        "simulate",
        help="Simulate and Label data ready for processing with the implemented models, e.g. for testing",
    )
    # unet_subparsers = unet_parsers.add_subparsers(dest="unet_subparsers")

    # Arguments for the general simulate command
    # parser = unet_subparsers.add_parser("simulate", help="simulate data for training")
    parser.add_argument(
        "--demes",
        type=existed_file,
        required=True,
        help="demographic model in the DEMES format",
    )
    parser.add_argument(
        "--nref",
        type=positive_int,
        required=True,
        help="number of samples in the reference population",
    )
    parser.add_argument(
        "--ntgt",
        type=positive_int,
        required=True,
        help="number of samples in the target population",
    )
    parser.add_argument(
        "--ref-id",
        type=str,
        required=True,
        help="name of the reference population in the demographic model",
        dest="ref_id",
    )
    parser.add_argument(
        "--tgt-id",
        type=str,
        required=True,
        help="name of the target population in the demographic model",
        dest="tgt_id",
    )
    parser.add_argument(
        "--src-id",
        type=str,
        required=True,
        help="name of the source population in the demographic model",
        dest="src_id",
    )
    parser.add_argument(
        "--seq-len",
        type=positive_int,
        required=True,
        help="length of the simulated genomes",
        dest="seq_len",
    )
    parser.add_argument(
        "--ploidy",
        type=positive_int,
        default=2,
        help="ploidy of the simulated genomes; default: 2",
    )
    parser.add_argument(
        "--phased",
        action="store_true",
        help="enable to use phased genotypes; default: False",
    )
    parser.add_argument(
        "--mut-rate",
        type=positive_number,
        default=1e-8,
        help="mutation rate per base pair per generation for the simulation; default: 1e-8",
        dest="mut_rate",
    )
    parser.add_argument(
        "--rec-rate",
        type=positive_number,
        default=1e-8,
        help="recombination rate per base pair per generation for the simulation; default: 1e-8",
        dest="rec_rate",
    )

    parser.add_argument(
        "--replicate",
        type=positive_int,
        default=1,
        help="total number of simulations, i.e. how many vcf-files with simulations are created; default: 1",
    )

    parser.add_argument(
        "--output-prefix",
        type=str,
        required=True,
        help="prefix of the output file name",
        dest="output_prefix",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="directory of the output files",
        dest="output_dir",
    )

    parser.add_argument(
        "--nprocess",
        type=positive_int,
        default=1,
        help="number of processesfor the simulation, i.e. how many vcf-files with simulations are processed in parallel; default: 1",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="random seed for the simulation; default: None",
    )

    # mutation model

    parser.add_argument(
        "--mutation_model",
        type=str,
        default="JC69",
        required=False,
        help="Mutation model to be used to for msprime simulations, currently binary ('binary') and Jukes-Cantor ('JC69', default) are directly supported",
        dest="mutation_model",
    )

    # mutation and recombination map parameters

    parser.add_argument(
        "--mut_map",
        type=str,
        default=None,
        required=False,
        help="either path to the mutation map file or instruction how to build it, currently 'gamma' and 'linear' are supported",
        dest="mut_map",
    )

    parser.add_argument(
        "--rec_map",
        type=str,
        default=None,
        required=False,
        help="either path to the recombination map file or instruction how to build it, currently 'linear' is supported",
        dest="rec_map",
    )

    parser.add_argument(
        "--batch_processing",
        action="store_true",
        help="whether to use the tree batch processing for in the function for obtaining true tracts",
        dest="true_tracts_batch_processing",
    )


    parser.add_argument(
        "--L_mut_map",
        default=100000,
        help="mean block length for mutation map construction (if chosen); default: 100000",
        dest="L_mut_map",
    )
    parser.add_argument(
        "--cv_mut_map",
        default=0.15,
        help="coefficient of variation for mutation map construction using exponential method (if chosen); default: 0.15",
        dest="cv_mut_map",
    )


    parser.add_argument(
        "--shuffle_rec_map",
        action="store_true",
        help="whether the rates of a loaded recombination rate map should also be shuffled",
        dest="shuffle_rec_map",
    )


    parser.add_argument(
        "--start_rep",
        type=int,
        default=0,
        help="first nr of the subfolders with replicates; default: 0",
    )

    parser.set_defaults(runner=_run_simulation)
