import demes
import msprime
import os
import tskit
import pyranges as pr
import math
from copy import deepcopy
from typing import Tuple
import pandas as pd
import pybedtools
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor as Pool


# version with batch processing


# functions
def create_final_tracts(tract_file):
    if tract_file is None or len(tract_file) == 0:
        # Return an empty DataFrame with expected schema
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "Sample"])

    tract_file.columns = ["Chromosome", "Start", "End", "Sample"]
    pr_tracts = pr.from_dict(tract_file)
    pr_tracts = pr_tracts.merge(strand=False, by="Sample")
    return pr_tracts


def partition_into_partitions(lst, num_partitions):
    partition_size = math.ceil(len(lst) / num_partitions)
    # partitions = [lst[i:i+partition_size] for i in range(0, len(lst), partition_size)]
    return [lst[i : i + partition_size] for i in range(0, len(lst), partition_size)]


# -----------------------------------------------------------------------------------------------------------------------

# def tree_batch_processing(args):
#    tree_batches, tree_object, migration_tree_object, bed_list, tgt_id = args
#    ts = tree_object.value
#    migration_ts = migration_tree_object.value
#    for m in migration_ts.migrations():
#        # Tree-sequences are sorted by the left ends of the intervals
#        # Can skip those tree-sequences are not overlapped with the interval of i.
#        for tree_index in tree_batches:
#            t = ts.at_index(tree_index)
#            if m.left >= t.interval.right: continue
#            if m.right <= t.interval.left: break # [l, r)
#            for n in ts.samples(tgt_id):
#                if t.is_descendant(n, m.node):
#                    left = m.left if m.left > t.interval.left else t.interval.left
#                    right = m.right if m.right < t.interval.right else t.interval.right
#                    bed_list.append([1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}'])
#    return True


# def migration_batch_processing(args):
#    tree_object, migration_batches, migration_tree_object, bed_list, tgt_id = args
#    ts = tree_object.value
#    migration_ts = migration_tree_object.value
#    for migration_index in migration_batches:
#        m = migration_ts.migrations()[migration_index]
#        # Tree-sequences are sorted by the left ends of the intervals
#        # Can skip those tree-sequences are not overlapped with the interval of i.
#        for t in ts.trees():
#            if m.left >= t.interval.right: continue
#            if m.right <= t.interval.left: break # [l, r)
#            for n in ts.samples(tgt_id):
#                if t.is_descendant(n, m.node):
#                    left = m.left if m.left > t.interval.left else t.interval.left
#                    right = m.right if m.right < t.interval.right else t.interval.right
#                    bed_list.append([1, int(left), int(right), f'tsk_{ts.node(n).individual}_{int(n%ploidy+1)}'])
#    return True
# -----------------------------------------------------------------------------------------------------------------------
# refactor tree_batch_processing & migration_batch_processing
def process_interval(ts, t, m, tgt_id, bed_list, ploidy, is_phased):
    if m.left >= t.interval.right:
        return
    if m.right <= t.interval.left:
        return  # exit early; caller should break if applicable
    for n in ts.samples(tgt_id):
        if t.is_descendant(n, m.node):
            left = m.left if m.left > t.interval.left else t.interval.left
            right = m.right if m.right < t.interval.right else t.interval.right
            if is_phased:
                bed_list.append(
                    [
                        1,
                        int(left),
                        int(right),
                        f"tsk_{ts.node(n).individual}_{int(n % ploidy + 1)}",
                    ]
                )
            else:
                bed_list.append(
                    [1, int(left), int(right), f"tsk_{ts.node(n).individual}"]
                )


def tree_batch_processing(args):

    (
        tree_batches,
        tree_object,
        migration_tree_object,
        bed_list,
        tgt_id,
        ploidy,
        is_phased,
    ) = args
    ts = tree_object.value
    migration_ts = migration_tree_object.value
    for m in migration_ts.migrations():
        for tree_index in tree_batches:
            t = ts.at_index(tree_index)
            if m.left >= t.interval.right:
                continue
            if m.right <= t.interval.left:
                break  # stop processing further trees for this migration
            process_interval(ts, t, m, tgt_id, bed_list, ploidy, is_phased)
    return True


def migration_batch_processing(args):

    (
        tree_object,
        migration_batches,
        migration_tree_object,
        bed_list,
        tgt_id,
        ploidy,
        is_phased,
    ) = args
    ts = tree_object.value
    migration_ts = migration_tree_object.value
    for migration_index in migration_batches:
        m = migration_ts.migrations()[migration_index]
        for t in ts.trees():
            if m.left >= t.interval.right:
                continue
            if m.right <= t.interval.left:
                break
            process_interval(ts, t, m, tgt_id, bed_list, ploidy, is_phased)
    return True


def get_true_tracts_batch(
    ts: tskit.TreeSequence,
    tgt_id: str,
    src_id: str,
    ploidy: int,
    is_phased: bool,
    bed_file=None,
) -> str:

    ### perhaps necessary!
    # mp.set_start_method("spawn", force=True)

    # index of population (source, target)
    src_name = [p.id for p in ts.populations() if p.metadata["name"] == src_id][0]
    tgt_name = [p.id for p in ts.populations() if p.metadata["name"] == tgt_id][0]

    # -----------------------------------------------------------------------------------------------------------------------
    # 1) extract ts with only relevant migrations
    ts_dump_mig = (
        ts.dump_tables()
    )  # make a copy of the tree sequence tables, for editing
    migtable = ts_dump_mig.migrations
    migtable2 = deepcopy(migtable)
    migtable2.clear()  # empty
    # we search for all rows involving source and target
    for mrow in migtable:
        if (mrow.dest == src_name) and (mrow.source == tgt_name):
            migtable2.append(mrow)
    # -----------------------------------------------------------------------------------------------------------------------
    # if no migrations b/n source & target, eg when simulating with no_introgression demes model to test the tools, could exit here
    # if migtable2.num_rows == 0:
    #    print("no migrations between source and target")
    #    sys.exit()
    ### commented: we want to create a truth tract file nonetheless...

    # -----------------------------------------------------------------------------------------------------------------------
    # the new tree sequence stores only the relevant migrations (involving source and target)
    ts_dump_mig.migrations.replace_with(migtable2)
    ts_dump_sequence_mig = ts_dump_mig.tree_sequence()
    # -----------------------------------------------------------------------------------------------------------------------

    # 2) remove all migrations & simplify tree
    ts_dump = ts.dump_tables()
    ts_dump.migrations.clear()
    ts_dump_sequence = ts_dump.tree_sequence()

    # we search for all nodes involving the relevant populations
    populations_not_to_remove = [src_name, tgt_name]
    # individuals_not_to_remove = []
    # for ind in ts.nodes():
    #    if ind.population in populations_not_to_remove:
    #        individuals_not_to_remove.append(ind.id)
    individuals_not_to_remove = [
        ind.id for ind in ts.nodes() if ind.population in populations_not_to_remove
    ]

    # the tree sequence object without migrations can be simplified
    # the simplification contains only the relevant (source-target involving) information
    ts_dump_sequence_simplified = ts_dump_sequence.simplify(
        individuals_not_to_remove,
        filter_populations=False,
        filter_individuals=False,
        filter_sites=False,
        filter_nodes=False,
    )

    # check has got here
    # print(ts_dump_sequence_simplified)

    tracts = pd.DataFrame(columns=["chrom", "start", "end", "sample"])
    tracts.columns = ["Chromosome", "Start", "End", "Sample"]

    manager = mp.Manager()  # create manager
    bed_list_tree_batches = (
        mp.Manager().list()
    )  # manager to write shared file - should this be inside or outside the run_pool function?

    # -----------------------------------------------------------------------------------------------------------------------
    # run pooling as a function, so can use either trees or migration fn
    def run_pool(x_inds, batch_function):
        pool = mp.Pool()  # initialise pool
        if len(x_inds) > 0:
            ## either partition across trees or across branches
            x_inds_batches = partition_into_partitions(x_inds, mp.cpu_count())
            ## make args
            args = [
                (
                    manager_trees,
                    curr_index,
                    manager_migrations,
                    bed_list_tree_batches,
                    tgt_name,
                    ploidy,
                    is_phased,
                )
                for curr_index in x_inds_batches
            ]
            all_done = pool.map(batch_function, args)
            all_done = list(all_done)
            # for concurrent.futures close and join do not exist, so in this case one would have to comment the next two lines
            pool.close()
            pool.join()
        return True

    # -----------------------------------------------------------------------------------------------------------------------
    migration_inds = list(range(ts_dump_sequence_mig.num_migrations))
    manager_migrations = manager.Value(tskit.trees.TreeSequence, ts_dump_sequence_mig)
    tree_inds = list(range(ts_dump_sequence_simplified.num_trees))
    manager_trees = manager.Value(tskit.trees.TreeSequence, ts_dump_sequence_simplified)

    # calc ratio between number of relevant migrations & number of trees
    # cutoff can be adjusted (0.5 is arbitrary choice atm)
    if (
        ts_dump_sequence_mig.num_migrations / ts_dump_sequence_simplified.num_trees
    ) < 0.5:
        # iterate over migrations, run_pool function
        all_done = run_pool(migration_inds, migration_batch_processing)
    else:
        # iterate over trees, run_pool function
        all_done = run_pool(tree_inds, tree_batch_processing)
    # -----------------------------------------------------------------------------------------------------------------------

    # now we have to convert the shared list into a pd dataframe
    list_bed_list_tree_batches = list(bed_list_tree_batches)
    bed_list_tree_df = pd.DataFrame(list_bed_list_tree_batches)
    bed_list_tree_df = create_final_tracts(bed_list_tree_df)
    try:
        bed_list_tree_df = bed_list_tree_df.df
    except:
        pass

    # check has got here
    # print(bed_list_tree_df)

    # writing as usual
    # print(f"bed_file: {bed_file}")
    if bed_file:
        open(bed_file, "w").close()

        if len(bed_list_tree_df) > 0:
            for s in bed_list_tree_df["Sample"].unique():
                # in case the simplified tree is used, there could be individuals set to -1
                if "-" not in s:
                    sample_tracts = bed_list_tree_df[bed_list_tree_df["Sample"] == s]
                    sample_tracts = (
                        pybedtools.BedTool.from_dataframe(sample_tracts)
                        .sort()
                        .merge()
                        .to_dataframe()
                    )
                    sample_tracts["Sample"] = s
                    sample_tracts.to_csv(
                        bed_file, sep="\t", mode="a", header=False, index=False
                    )
                else:
                    # print("no migration events found due to wrong individual id")
                    df = pd.DataFrame()
                    df.to_csv(bed_file, sep="\t", mode="a", header=False, index=False)

        else:
            print("no migration events found, the bed list is empty")
            df = pd.DataFrame()

            df.to_csv(bed_file, sep="\t", mode="a", header=False, index=False)

        return True
    return bed_list_tree_df