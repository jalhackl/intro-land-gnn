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


import os
import random
import shutil
import pandas as pd
from multiprocessing import Lock, Value
from gaishisim.utils.multiprocessing import mp_manager
from gaishisim.utils.generators import RandomNumberGenerator
from gaishisim.utils.simulators import MsprimeSimulator



def simulate_test_data(
    demo_model_file: str,
    nrep: int,
    nref: int,
    ntgt: int,
    ref_id: str,
    tgt_id: str,
    src_id: str,
    ploidy: int,
    is_phased: bool,
    seq_len: int,
    mut_rate: float,
    rec_rate: float,
    output_prefix: str,
    output_dir: str,
    seed: int,
    nprocess: int,
    mutation_model: str,
    rec_map: str,
    mut_map: str,
    L_mut_map = 100000,
    cv_mut_map = 0.15,
    shuffle_rec_map: bool = False,
    true_tracts_batch_processing: bool = False,
    start_rep=0
) -> None:
    """
    Simulates genomic data using a demographic model and processes it for downstream analysis.

    This function utilizes `MsprimeSimulator` to generate synthetic genomic data and
    `RandomNumberGenerator` to control the simulation process. It then employs a
    multiprocessing manager to efficiently process the generated data.

    Args:
        demo_model_file (str): Path to the demographic model file used for simulations.
        nrep (int): Number of simulation replicates to generate.
        nref (int): Number of reference individuals included in the simulation.
        ntgt (int): Number of target individuals included in the simulation.
        ref_id (str): Identifier for the reference population in the simulation.
        tgt_id (str): Identifier for the target population in the simulation.
        src_id (str): Identifier for the source population in the simulation.
        ploidy (int): Ploidy level of the individuals in the simulation (e.g., 2 for diploids).
        is_phased (bool): If True, the simulated data is phased; otherwise, it is unphased.
        seq_len (int): Length of the genomic sequence to be simulated (in base pairs).
        mut_rate (float): Mutation rate per base pair per generation.
        rec_rate (float): Recombination rate per base pair per generation.
        output_prefix (str): Prefix for output file names.
        output_dir (str): Directory where output files will be stored.
        seed (int): Random seed for reproducibility.
        nprocess (int): Number of processes for parallel processing.
        mutation_model (str): Mutation model to be used to for msprime simulations, currently binary and Jukes-Cantor are directly supported.
        rec_map (str)
        mut_map (str)
        shuffle_rec_map (str)
        true_tracts_batch_processing (bool)
        
    Returns:
        None: The function does not return anything but generates simulation output files.

    Raises:
        ValueError: If an invalid parameter is provided (e.g., negative mutation rate).
        SystemExit: If an error occurs within the multiprocessing manager.

    Notes:
        - The function integrates with `MsprimeSimulator` for data generation.
        - It ensures reproducibility using a controlled random seed.
        - Outputs are saved in the specified directory.
        - Parallel processing is leveraged to speed up simulations.
    """

    # Initialize the simulator with the given parameters

    simulator = MsprimeSimulator(
        demo_model_file=demo_model_file,
        nref=nref,
        ntgt=ntgt,
        ref_id=ref_id,
        tgt_id=tgt_id,
        src_id=src_id,
        ploidy=ploidy,
        seq_len=seq_len,
        mut_rate=mut_rate,
        rec_rate=rec_rate,
        output_prefix=output_prefix,
        output_dir=output_dir,
        is_phased=is_phased,
        mutation_model=mutation_model,
        rec_map=rec_map,
        mut_map=mut_map,
        L_mut_map = L_mut_map,
        cv_mut_map = cv_mut_map,
        shuffle_rec_map=shuffle_rec_map,
        true_tracts_batch_processing=true_tracts_batch_processing,
    )

    generator = RandomNumberGenerator(
        start_rep=start_rep,
        nrep=nrep,
        seed=seed,
    )

    mp_manager(
        job=simulator,
        data_generator=generator,
        nprocess=nprocess,
    )
