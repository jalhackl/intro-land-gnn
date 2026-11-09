"""
Modified version of `convert_to_tskit.py` from the SINGER github repository
(github.com/popgenmethods/SINGER) that computes mutation parents correctly so
as to load without error in tskit>=1.0
"""

import sys
import os
import argparse
import numpy as np
import tskit


def read_ts(node_file, edge_file):
    node_time = np.loadtxt(node_file)
    edge_span = np.loadtxt(edge_file)
    edge_span = edge_span[edge_span[:, 2] >= 0, :]
    length = max(edge_span[:, 1])
    tables = tskit.TableCollection(sequence_length=length)
    node_table = tables.nodes
    edge_table = tables.edges
    prev_time = -1
    for t in node_time:
        if (t == 0):
            node_table.add_row(flags=tskit.NODE_IS_SAMPLE)
        else:
            t = max(prev_time + 1e-4, t)
            node_table.add_row(time = t)
            prev_time = t
    parent_indices = np.array(edge_span[:, 2], dtype = np.int32)
    child_indices = np.array(edge_span[:, 3], dtype = np.int32)
    edge_table.set_columns(left = edge_span[:, 0], right = edge_span[:, 1], parent = parent_indices, child = child_indices)
    return tables


def read_mutation(tables, mutation_file, apply_zero_convention=False, remove_invalid_mutations=True):
    mutations = np.loadtxt(mutation_file)

    # Remove mutations outside the valid sequence interval
    if not apply_zero_convention:
        invalid = (
            (mutations[:, 0] < 0) |
            (mutations[:, 0] >= tables.sequence_length)
        )
    else:
        invalid = (
            (mutations[:, 0] < 1) |
            (mutations[:, 0] > tables.sequence_length)
        )

    if np.any(invalid):
        print(
            f"WARNING: {np.sum(invalid)} mutations "
            f"outside sequence boundaries."
        )
        print("Sequence length:", tables.sequence_length)
        print("Invalid mutations:", mutations[invalid])


    if remove_invalid_mutations:
        print("Removing invalid mutations.")
        mutations = mutations[~invalid]


    n = mutations.shape[0]
    mut_pos = 0
    for i in range(n):
        if mutations[i, 0] != mut_pos:
            if not apply_zero_convention:
                tables.sites.add_row(position=mutations[i, 0], ancestral_state='0')
            else:
                tables.sites.add_row(position=mutations[i, 0]-1, ancestral_state='0')
            mut_pos = mutations[i, 0]
        site_id = tables.sites.num_rows - 1
        tables.mutations.add_row(site=site_id, node=int(mutations[i, 1]), derived_state=str(int(mutations[i, 3])))
    # added: rebuild mutations table in reverse time order at each position
    mut_time = tables.nodes.time[tables.mutations.node]
    mut_coord = tables.sites.position[tables.mutations.site]
    mut_order = np.lexsort((-mut_time, mut_coord))
    mut_state = tskit.unpack_strings(
        tables.mutations.derived_state, 
        tables.mutations.derived_state_offset,
    )
    mut_state, mut_state_offset = tskit.pack_strings(np.array(mut_state)[mut_order])
    tables.mutations.set_columns(
        site=tables.mutations.site[mut_order],
        node=tables.mutations.node[mut_order],
        time=np.repeat(tskit.UNKNOWN_TIME, tables.mutations.num_rows),
        derived_state=mut_state,
        derived_state_offset=mut_state_offset,
    )


def read_ARG(node_file, branch_file, mutation_file, apply_zero_convention=False, remove_invalid_mutations=True):
    tables = read_ts(node_file, branch_file)
    read_mutation(tables, mutation_file, apply_zero_convention=apply_zero_convention, remove_invalid_mutations=remove_invalid_mutations)
    # added: calculate mutation parents and times
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    ts = tables.tree_sequence()
    return ts


def write_fast_trees(input_prefix, output_prefix, start, end, step, apply_zero_convention=False, remove_invalid_mutations=True):
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}_fast_nodes_{i}.txt"
        branch_file = f"{input_prefix}_fast_branches_{i}.txt"
        mutation_file = f"{input_prefix}_fast_muts_{i}.txt"
        ts = read_ARG(node_file, branch_file, mutation_file, apply_zero_convention=apply_zero_convention, remove_invalid_mutations=remove_invalid_mutations)
        ts.dump(trees_file)


def write_trees(input_prefix, output_prefix, start, end, step, apply_zero_convention=False, remove_invalid_mutations=True):
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}_nodes_{i}.txt"
        branch_file = f"{input_prefix}_branches_{i}.txt"
        mutation_file = f"{input_prefix}_muts_{i}.txt"
        ts = read_ARG(node_file, branch_file, mutation_file, apply_zero_convention=apply_zero_convention, remove_invalid_mutations=remove_invalid_mutations)
        ts.dump(trees_file) 

