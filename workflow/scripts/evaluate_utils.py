import pybedtools
import numpy as np
import os



import pandas as pd
import pyranges as pr
import numpy as np
import os
#from scipy.special import expit
#from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
#import pickle
#import h5py
from typing import Optional, Tuple, List
from copy import deepcopy

# evaluation function from gaia beta, currently not used
def window_evaluate(
    true_tract_file: str,
    inferred_tract_file: str,
    seq_len: int,
    sample_size: int,
    ploidy: int,
    is_phased: bool,
    output: str,
    column_names = None
) -> None:
    """
    Evaluates model performance with precision and recall based on genomic windows by comparing
    true and inferred introgressed tracts.

    Parameters
    ----------
    true_tract_file : str
        File path for the BED file containing true introgressed fragments.
    inferred_tract_file : str
        File path for the BED file containing inferred introgressed fragments.
    seq_len : int
        The total length of the sequence in base pairs.
    sample_size : int
        The number of samples analyzed in the study.
    ploidy : int
        The ploidy of the genomes being analyzed.
    is_phased : bool
        Indicates whether the genetic data is phased.
    output : str
        File path for the output file storing the model performance metrics in a tab-separated format.

    """
    if not column_names:
        column_names = ["Chromosome", "Start", "End", "Sample"]

    chrom_column = column_names[0]
    start_column = column_names[1]
    end_column = column_names[2]
    sample_column = column_names[3]


    try:
        true_tracts = pd.read_csv(
            true_tract_file,
            sep="\t",
            header=None,
            names=column_names,
        )
    except pd.errors.EmptyDataError:
        true_tracts_samples = []
    else:
        true_tracts_samples = true_tracts[sample_column].unique()
        true_tracts = pr.PyRanges(true_tracts).merge(by=sample_column)

    try:
        inferred_tracts = pd.read_csv(
            inferred_tract_file,
            sep="\t",
            header=None,
            names=column_names,
        )
        inferred_tracts[end_column] = inferred_tracts[end_column].clip(upper=seq_len)
    except pd.errors.EmptyDataError:
        inferred_tracts_samples = []
    else:
        inferred_tracts_samples = inferred_tracts[sample_column].unique()
        inferred_tracts = pr.PyRanges(inferred_tracts).merge(by=sample_column)

    if is_phased:
        sample_size = sample_size * ploidy

    res = pd.DataFrame(
        columns=[
            "Sample",
            "Sequence_length",
            "True_tracts_length",
            "Inferred_tracts_length",
            "True_positives_length",
            "False_positives_length",
            "True_negatives_length",
            "False_negatives_length",
        ]
    )

    # sum_ntrue_tracts = 0
    # sum_ninferred_tracts = 0
    # sum_ntrue_positives = 0

    for s in np.intersect1d(true_tracts_samples, inferred_tracts_samples):
        ind_true_tracts = true_tracts[true_tracts.Sample == s]
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]
        ind_overlaps = ind_true_tracts.intersect(ind_inferred_tracts)

        ntrue_tracts = np.sum(
            [x[1].End.astype("int") - x[1].Start.astype("int") for x in ind_true_tracts]
        )
        ninferred_tracts = np.sum(
            [
                x[1].End.astype("int") - x[1].Start.astype("int")
                for x in ind_inferred_tracts
            ]
        )
        ntrue_positives = np.sum(
            [x[1].End.astype("int") - x[1].Start.astype("int") for x in ind_overlaps]
        )
        nfalse_positives = ninferred_tracts - ntrue_positives
        nfalse_negatives = ntrue_tracts - ntrue_positives
        ntrue_negatives = seq_len - ntrue_tracts - nfalse_positives
        # precision, recall = cal_pr(ntrue_tracts, ninferred_tracts, ntrue_positives)
        res.loc[len(res.index)] = [
            s,
            seq_len,
            ntrue_tracts,
            ninferred_tracts,
            ntrue_positives,
            nfalse_positives,
            ntrue_negatives,
            nfalse_negatives,
        ]

        # sum_ntrue_tracts += ntrue_tracts
        # sum_ninferred_tracts += ninferred_tracts
        # sum_ntrue_positives += ntrue_positives

    for s in np.setdiff1d(true_tracts_samples, inferred_tracts_samples):
        # ninferred_tracts = 0
        ind_true_tracts = true_tracts[true_tracts.Sample == s]
        ntrue_tracts = np.sum(
            [x[1].End.astype("int") - x[1].Start.astype("int") for x in ind_true_tracts]
        )
        res.loc[len(res.index)] = [
            s,
            seq_len,
            ntrue_tracts,
            0,
            0,
            0,
            seq_len,
            ntrue_tracts,
        ]

        # sum_ntrue_tracts += ntrue_tracts

    for s in np.setdiff1d(inferred_tracts_samples, true_tracts_samples):
        # ntrue_tracts = 0
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]
        ninferred_tracts = np.sum(
            [
                x[1].End.astype("int") - x[1].Start.astype("int")
                for x in ind_inferred_tracts
            ]
        )
        res.loc[len(res.index)] = [
            s,
            seq_len,
            0,
            ninferred_tracts,
            0,
            ninferred_tracts,
            seq_len - ninferred_tracts,
            0,
        ]

        # sum_ninferred_tracts += ninferred_tracts

    # sum_nfalse_positives = sum_ninferred_tracts - sum_ntrue_positives
    # sum_nfalse_negatives = sum_ntrue_tracts - sum_ntrue_positives
    # sum_ntrue_negatives = seq_len * sample_size - sum_ntrue_tracts - sum_nfalse_positives

    res = res.sort_values(by=["Sample"])

    numeric_columns = res.select_dtypes(include=[float, int]).columns
    column_means = res[numeric_columns].mean()
    mean_df = pd.DataFrame([column_means], index=["Mean"])
    mean_df.insert(0, "Sample", "Average")
    res = pd.concat([res, mean_df])

    # total_precision, total_recall = cal_pr(sum_ntrue_tracts, sum_ninferred_tracts, sum_ntrue_positives)
    # total_len = seq_len*sample_size
    # res.loc[len(res.index)] = [
    #    'Total',
    #    total_len,
    #    sum_ntrue_tracts,
    #    sum_ninferred_tracts,
    #    sum_ntrue_positives/sum_ntrue_tracts*100 if sum_ntrue_tracts != 0 else 0,
    #    sum_nfalse_positives/(total_len-sum_ntrue_tracts)*100,
    #    sum_ntrue_negatives/(total_len-sum_ntrue_tracts)*100,
    #    sum_nfalse_negatives/sum_ntrue_tracts*100 if sum_ntrue_tracts != 0 else 0,
    # ]

    res.fillna("NaN").to_csv(output, sep="\t", index=False)



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

    try:
        truth_tracts = pybedtools.BedTool(true_tracts).sort().merge()
        inferred_tracts =  pybedtools.BedTool(inferred_tracts).sort().merge()

        total_inferred_tracts = sum(x.stop - x.start for x in (inferred_tracts))
        total_true_tracts =  sum(x.stop - x.start for x in (truth_tracts))
        true_positives = sum(x.stop - x.start for x in inferred_tracts.intersect(truth_tracts))

        if float(total_inferred_tracts) == 0: precision = np.nan
        else: precision = true_positives / float(total_inferred_tracts) * 100
        if float(total_true_tracts) == 0: recall = np.nan
        else: recall = true_positives / float(total_true_tracts) * 100

        return precision, recall

    except Exception as e:
        return 0, 0


def cal_accuracy_samples(true_tracts, inferred_tracts, column_names=None):
    """
    Calculate base-pair precision and recall for a single pair of BED files,
    respecting a 'sample' column.

    Expected BED columns:
        chrom, start, end, sample
    """

    import pandas as pd
    import numpy as np
    import pybedtools

    if not column_names:
        column_names = ["chrom", "start", "end", "sample"]
    chrom_column = column_names[0]
    start_column = column_names[1]
    end_column = column_names[2]
    sample_column = column_names[3]


    # Read BED files
    truth_df = pd.read_csv(
        true_tracts, sep="\t", header=None,
        names=column_names
    )

    inferred_df = pd.read_csv(
        inferred_tracts, sep="\t", header=None,
        names=column_names
    )

    total_true = 0
    total_inferred = 0
    total_true_positive = 0

    # Get all samples present in either file
    samples = set(truth_df[sample_column]).union(set(inferred_df[sample_column]))

    for s in samples:

        truth_s = truth_df[truth_df[sample_column] == s][[chrom_column, start_column, end_column]]
        inferred_s = inferred_df[inferred_df[sample_column] == s][[chrom_column, start_column, end_column]]

        # Convert to BedTool if non-empty
        if not truth_s.empty:
            truth_bt = pybedtools.BedTool.from_dataframe(truth_s).sort().merge()
            t_true = sum(x.stop - x.start for x in truth_bt)
        else:
            t_true = 0

        if not inferred_s.empty:
            inferred_bt = pybedtools.BedTool.from_dataframe(inferred_s).sort().merge()
            t_inf = sum(x.stop - x.start for x in inferred_bt)
        else:
            t_inf = 0

        if t_true > 0 and t_inf > 0:
            t_tp = sum(x.stop - x.start
                       for x in inferred_bt.intersect(truth_bt))
        else:
            t_tp = 0

        total_true += t_true
        total_inferred += t_inf
        total_true_positive += t_tp

    precision = (
        total_true_positive / total_inferred * 100
        if total_inferred > 0 else np.nan
    )

    recall = (
        total_true_positive / total_true * 100
        if total_true > 0 else np.nan
    )

    return precision, recall


def cal_accuracy_multiple_samples(true_tract_files, inferred_tract_files, column_names=None):
    """
    Calculate aggregated base-pair precision and recall across
    multiple BED file pairs with a 'sample' column.

    Expected BED columns:
        chrom, start, end, sample
    """

    import pandas as pd
    import numpy as np
    import pybedtools
    import os

    if not column_names:
        column_names = ["chrom", "start", "end", "sample"]
    chrom_column = column_names[0]
    start_column = column_names[1]
    end_column = column_names[2]
    sample_column = column_names[3]



    total_true = 0
    total_inferred = 0
    total_true_positive = 0

    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):

        # Handle missing or empty files safely
        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth_df = pd.DataFrame(columns=column_names)
        else:
            truth_df = pd.read_csv(
                true_f, sep="\t", header=None,
                names=column_names
            )

        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred_df = pd.DataFrame(columns=column_names)
        else:
            inferred_df = pd.read_csv(
                inf_f, sep="\t", header=None,
                names=column_names
            )

        samples = set(truth_df[sample_column]).union(set(inferred_df[sample_column]))

        for s in samples:

            truth_s = truth_df[truth_df[sample_column] == s][[chrom_column, start_column, end_column]]
            inferred_s = inferred_df[inferred_df[sample_column] == s][[chrom_column, start_column, end_column]]

            if not truth_s.empty:
                truth_bt = pybedtools.BedTool.from_dataframe(truth_s).sort().merge()
                t_true = sum(x.stop - x.start for x in truth_bt)
            else:
                t_true = 0

            if not inferred_s.empty:
                inferred_bt = pybedtools.BedTool.from_dataframe(inferred_s).sort().merge()
                t_inf = sum(x.stop - x.start for x in inferred_bt)
            else:
                t_inf = 0

            if t_true > 0 and t_inf > 0:
                t_tp = sum(x.stop - x.start
                           for x in inferred_bt.intersect(truth_bt))
            else:
                t_tp = 0

            total_true += t_true
            total_inferred += t_inf
            total_true_positive += t_tp

    precision = (
        total_true_positive / total_inferred * 100
        if total_inferred > 0 else np.nan
    )

    recall = (
        total_true_positive / total_true * 100
        if total_true > 0 else np.nan
    )

    return precision, recall


def cal_accuracy_multiple(true_tract_files, inferred_tract_files):
    """
    Calculate precision and recall across many tract files
    by aggregating base-pair lengths before computing metrics.

    Arguments:
        true_tract_files      list[str]: BED files containing true introgressed tracts.
        inferred_tract_files  list[str]: BED files containing inferred introgressed tracts.

    Returns:
        precision float
        recall float
    """
    total_true = 0
    total_inferred = 0
    total_true_positive = 0

    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):

        # Skip missing or empty inputs safely
        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth = pybedtools.BedTool()  # empty
            t_true = 0
        else:
            truth = pybedtools.BedTool(true_f).sort().merge()
            t_true = sum(x.stop - x.start for x in truth)


        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred = pybedtools.BedTool()
            t_inf = 0
        else:
            inferred = pybedtools.BedTool(inf_f).sort().merge()
            t_inf = sum(x.stop - x.start for x in inferred)

        # Compute lengths
        if t_inf == 0 or t_true == 0:
            t_tp = 0
        else:
            t_tp  = sum(x.stop - x.start for x in inferred.intersect(truth))

        total_true += t_true
        total_inferred += t_inf
        total_true_positive += t_tp

    # Compute aggregated precision and recall
    precision = (total_true_positive / total_inferred * 100
                 if total_inferred > 0 else np.nan)

    recall = (total_true_positive / total_true * 100
              if total_true > 0 else np.nan)

    return precision, recall




def cal_confusion_multiple(true_tract_files, inferred_tract_files):
    """
    Calculate base-pair confusion matrix components (TP, FP, FN)
    across many tract files by aggregating lengths.

    Arguments:
        true_tract_files      list[str]: BED files with true introgressed tracts.
        inferred_tract_files  list[str]: BED files with inferred introgressed tracts.

    Returns:
        total_true_positive (TP)
        total_false_positive (FP)
        total_false_negative (FN)
    """

    total_true_positive = 0
    total_false_positive = 0
    total_false_negative = 0

    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):

        # Load truth
        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth = pybedtools.BedTool()
        else:
            truth = pybedtools.BedTool(true_f).sort().merge()

        # Load inferred
        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred = pybedtools.BedTool()
        else:
            inferred = pybedtools.BedTool(inf_f).sort().merge()

        # --- TP ---
        inter = inferred.intersect(truth).sort().merge()
        t_tp = sum(x.stop - x.start for x in inter)

        # --- FP = inferred minus truth ---
        fp = inferred.subtract(truth).sort().merge()
        t_fp = sum(x.stop - x.start for x in fp)

        # --- FN = truth minus inferred ---
        fn = truth.subtract(inferred).sort().merge()
        t_fn = sum(x.stop - x.start for x in fn)

        total_true_positive += t_tp
        total_false_positive += t_fp
        total_false_negative += t_fn

    return total_true_positive, total_false_positive, total_false_negative




def cal_metrics_multiple_wo_length_check(true_tract_files, inferred_tract_files, sequence_length=None, sequence_length_per_file=True):
    """
    Calculate precision, recall, and base-pair confusion matrix components
    across many tract files by aggregating base-pair lengths.

    Arguments:
        true_tract_files      list[str]: BED files containing true introgressed tracts.
        inferred_tract_files  list[str]: BED files containing inferred introgressed tracts.
        sequence_length       int or None: total length of the sequence/genome/chromosome.
                                          If provided, TN will be calculated; otherwise TN=None.
        sequence_length_per_file    bool: if true, sequence length is length of one file (all files are assumed to have the same nr of bases)

    Returns:
        precision float
        recall float
        total_true_positive (TP)
        total_false_positive (FP)
        total_false_negative (FN)
        total_true_negative (TN or None)
    """

    total_true = 0
    total_inferred = 0
    total_true_positive = 0
    total_false_positive = 0
    total_false_negative = 0

    if sequence_length is not None:
        if sequence_length_per_file:
            total_bases = 0
        else:
            total_bases = sequence_length

    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):

        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth = pybedtools.BedTool()
            t_true = 0
        else:
            truth = pybedtools.BedTool(true_f).sort().merge()
            t_true = sum(x.stop - x.start for x in truth)

        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred = pybedtools.BedTool()
            t_inf = 0
        else:
            inferred = pybedtools.BedTool(inf_f).sort().merge()
            t_inf = sum(x.stop - x.start for x in inferred)

        # --- True Positive ---
        if t_inf == 0 or t_true == 0:
            t_tp = 0
        else:
            inter = inferred.intersect(truth).sort().merge()
            t_tp = sum(x.stop - x.start for x in inter)

        # --- False Positive (inferred but not true) ---
        if t_true == 0:
            t_fp = t_inf
        elif t_inf == 0:
            t_fp = 0
        else:
            fp = inferred.subtract(truth).sort().merge()
            t_fp = sum(x.stop - x.start for x in fp)

        # --- False Negative (true but not inferred) ---
        if t_true == 0:
            t_fn = 0
        elif t_inf == 0:
            t_fn = t_true
        else:
            fn = truth.subtract(inferred).sort().merge()
            t_fn = sum(x.stop - x.start for x in fn)

        # Aggregate across files
        total_true += t_true
        total_inferred += t_inf
        total_true_positive += t_tp
        total_false_positive += t_fp
        total_false_negative += t_fn
        if sequence_length is not None:
            if sequence_length_per_file:
                total_bases += sequence_length

    precision = (total_true_positive / total_inferred * 100
                 if total_inferred > 0 else np.nan)

    recall = (total_true_positive / total_true * 100
              if total_true > 0 else np.nan)

    # --- True Negative (optional) ---
    if sequence_length is not None:
        used = total_true_positive + total_false_positive + total_false_negative
        total_true_negative = total_bases - used

    else:
        total_true_negative = None

    return precision,recall,total_true_positive,total_false_positive,total_false_negative,total_true_negative



def cal_metrics_multiple_no_samples(true_tract_files,
                         inferred_tract_files,
                         sequence_length=None,
                         sequence_length_per_file=True,
                         return_total_for_processing=False):

    total_true = 0
    total_inferred = 0
    total_true_positive = 0
    total_false_positive = 0
    total_false_negative = 0

    if sequence_length is not None:
        if sequence_length_per_file:
            total_bases = 0
        else:
            total_bases = sequence_length

    has_true = True
    has_inferred = True
    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):

        # ----------------------------
        # Load TRUE tracts
        # ----------------------------
        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth = pybedtools.BedTool()
            t_true = 0
            has_true = False
        else:
            truth = pybedtools.BedTool(true_f).sort().merge()
            t_true = sum(x.stop - x.start for x in truth)
        # ----------------------------
        # Load INFERRED tracts
        # ----------------------------
        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred = pybedtools.BedTool()
            t_inf = 0
            has_inferred = False
        else:
            inferred = pybedtools.BedTool(inf_f).sort().merge()
            t_inf = sum(x.stop - x.start for x in inferred)
        # ----------------------------
        # Clip to sequence_length (if provided)
        # ----------------------------
        if sequence_length is not None:

            # Collect chromosome names present in either file (however, only one chrom should be in the file)
            chroms = set()
            if has_true:
                for interval in truth:
                    chroms.add(interval.chrom)

            if has_inferred:
                for interval in inferred:
                    chroms.add(interval.chrom)

            if chroms:
                # Create clipping BED dynamically
                clip_string = ""
                for chrom in chroms:
                    clip_string += f"{chrom}\t0\t{sequence_length}\n"

                clip_bed = pybedtools.BedTool(clip_string, from_string=True)

                if has_true:
                    truth = truth.intersect(clip_bed).sort().merge()
                if has_inferred:
                    inferred = inferred.intersect(clip_bed).sort().merge()

        # Recompute lengths after clipping
        if has_true:
            t_true = sum(x.stop - x.start for x in truth)
        if has_inferred:
            t_inf = sum(x.stop - x.start for x in inferred)

        # ----------------------------
        # True Positive
        # ----------------------------
        if t_inf == 0 or t_true == 0:
            t_tp = 0
        else:
            inter = inferred.intersect(truth).sort().merge()
            t_tp = sum(x.stop - x.start for x in inter)

        # ----------------------------
        # False Positive
        # ----------------------------
        if t_true == 0:
            t_fp = t_inf
        elif t_inf == 0:
            t_fp = 0
        else:
            fp = inferred.subtract(truth).sort().merge()
            t_fp = sum(x.stop - x.start for x in fp)

        # ----------------------------
        # False Negative
        # ----------------------------
        if t_true == 0:
            t_fn = 0
        elif t_inf == 0:
            t_fn = t_true
        else:
            fn = truth.subtract(inferred).sort().merge()
            t_fn = sum(x.stop - x.start for x in fn)

        # ----------------------------
        # Aggregate
        # ----------------------------
        total_true += t_true
        total_inferred += t_inf
        total_true_positive += t_tp
        total_false_positive += t_fp
        total_false_negative += t_fn

        if sequence_length is not None and sequence_length_per_file:
            total_bases += sequence_length

    # ----------------------------
    # Precision & Recall
    # ----------------------------
    precision = (total_true_positive / total_inferred * 100
                 if total_inferred > 0 else np.nan)

    recall = (total_true_positive / total_true * 100
              if total_true > 0 else np.nan)

    # ----------------------------
    # True Negative (optional)
    # ----------------------------
    if sequence_length is not None:
        used = (total_true_positive +
                total_false_positive +
                total_false_negative)
        total_true_negative = total_bases - used
    else:
        total_true_negative = None

    if return_total_for_processing:
        return total_true, total_inferred, total_true_positive, total_false_positive, total_false_negative, total_true_negative

    return (precision,
            recall,
            total_true_positive,
            total_false_positive,
            total_false_negative,
            total_true_negative)


def cal_metrics_multiple(true_tract_files,
                                      inferred_tract_files,
                                      sequence_length=None,
                                      sequence_length_per_file=True,
                                      return_total_for_processing=False,
                                      column_names=None):

    import pandas as pd
    import numpy as np
    import pybedtools
    import os

    print("TRUE FILES:", true_tract_files)
    print("INFERRED FILES:", inferred_tract_files)

    if not column_names:
        column_names = ["chrom", "start", "end", "sample"]

    chrom_column = column_names[0]
    start_column = column_names[1]
    end_column = column_names[2]
    sample_column = column_names[3]

    total_true = 0
    total_inferred = 0
    total_true_positive = 0
    total_false_positive = 0
    total_false_negative = 0

    if sequence_length is not None:
        if sequence_length_per_file:
            total_bases = 0
        else:
            total_bases = sequence_length

    for true_f, inf_f in zip(true_tract_files, inferred_tract_files):
        print("CurrTRUE FILES:", true_f)
        print("CurrINFERRED FILES:", inf_f)
        # --- Load safely ---
        if not os.path.exists(true_f) or os.path.getsize(true_f) == 0:
            truth_df = pd.DataFrame(columns=column_names)
        else:
            truth_df = pd.read_csv(true_f, sep="\t", header=None,
                                   names=column_names)

        if not os.path.exists(inf_f) or os.path.getsize(inf_f) == 0:
            inferred_df = pd.DataFrame(columns=column_names)
        else:
            inferred_df = pd.read_csv(inf_f, sep="\t", header=None,
                                      names=column_names)

        samples = set(truth_df[sample_column]).union(
                  set(inferred_df[sample_column]))

        for s in samples:

            truth_s = truth_df[truth_df[sample_column] == s][
                [chrom_column, start_column, end_column]
            ]

            inferred_s = inferred_df[inferred_df[sample_column] == s][
                [chrom_column, start_column, end_column]
            ]

            # --- Convert to BedTool and merge ---
            if not truth_s.empty:
                truth_bt = (
                    pybedtools.BedTool.from_dataframe(truth_s)
                    .sort()
                    .merge()
                )
                t_true = sum(x.stop - x.start for x in truth_bt)
            else:
                truth_bt = None
                t_true = 0

            if not inferred_s.empty:
                inferred_bt = (
                    pybedtools.BedTool.from_dataframe(inferred_s)
                    .sort()
                    .merge()
                )
                t_inf = sum(x.stop - x.start for x in inferred_bt)
            else:
                inferred_bt = None
                t_inf = 0

            # --- Clip if sequence_length provided ---
            if sequence_length is not None:

                chroms = set()

                if truth_bt:
                    for interval in truth_bt:
                        chroms.add(interval.chrom)

                if inferred_bt:
                    for interval in inferred_bt:
                        chroms.add(interval.chrom)

                if chroms:
                    clip_string = ""
                    for chrom in chroms:
                        clip_string += f"{chrom}\t0\t{sequence_length}\n"

                    clip_bed = pybedtools.BedTool(
                        clip_string, from_string=True
                    )

                    if truth_bt:
                        truth_bt = truth_bt.intersect(clip_bed).sort().merge()
                        t_true = sum(x.stop - x.start for x in truth_bt)

                    if inferred_bt:
                        inferred_bt = inferred_bt.intersect(clip_bed).sort().merge()
                        t_inf = sum(x.stop - x.start for x in inferred_bt)

            # --- True Positive ---
            if t_true > 0 and t_inf > 0:
                inter = inferred_bt.intersect(truth_bt).sort().merge()
                t_tp = sum(x.stop - x.start for x in inter)
            else:
                t_tp = 0

            # --- False Positive ---
            if t_true == 0:
                t_fp = t_inf
            elif t_inf == 0:
                t_fp = 0
            else:
                fp = inferred_bt.subtract(truth_bt).sort().merge()
                t_fp = sum(x.stop - x.start for x in fp)

            # --- False Negative ---
            if t_true == 0:
                t_fn = 0
            elif t_inf == 0:
                t_fn = t_true
            else:
                fn = truth_bt.subtract(inferred_bt).sort().merge()
                t_fn = sum(x.stop - x.start for x in fn)

            # --- Aggregate ---
            total_true += t_true
            total_inferred += t_inf
            total_true_positive += t_tp
            total_false_positive += t_fp
            total_false_negative += t_fn

            if sequence_length is not None and sequence_length_per_file:
                total_bases += sequence_length

    # --- Precision & Recall ---
    precision = (
        total_true_positive / total_inferred * 100
        if total_inferred > 0 else np.nan
    )

    recall = (
        total_true_positive / total_true * 100
        if total_true > 0 else np.nan
    )

    # --- True Negative ---
    if sequence_length is not None:
        used = (total_true_positive +
                total_false_positive +
                total_false_negative)
        total_true_negative = total_bases - used
    else:
        total_true_negative = None

    if return_total_for_processing:
        return (total_true,
                total_inferred,
                total_true_positive,
                total_false_positive,
                total_false_negative,
                total_true_negative)

    return (precision,
            recall,
            total_true_positive,
            total_false_positive,
            total_false_negative,
            total_true_negative)




def cal_metrics_single(
    true_tract_file,
    inferred_tract_file,
    sequence_length=None,
    column_names=None,
    return_tn=False,
    return_precision_recall=False
):
    import pandas as pd
    import numpy as np
    import pybedtools
    import os

    if column_names is None:
        column_names = ["chrom", "start", "end", "sample"]

    chrom, start, end, sample = column_names

    # --- Load files safely ---
    def load(path):
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return pd.DataFrame(columns=column_names)
        return pd.read_csv(path, sep="\t", header=None, names=column_names)

    truth_df = load(true_tract_file)
    inferred_df = load(inferred_tract_file)

    # --- get samples ---
    samples = set(truth_df[sample]).union(set(inferred_df[sample]))

    total_tp = total_fp = total_fn = 0
    total_true = total_inf = 0

    total_bases = sequence_length if sequence_length is not None else 0

    for s in samples:
        truth_s = truth_df[truth_df[sample] == s][[chrom, start, end]]
        inf_s = inferred_df[inferred_df[sample] == s][[chrom, start, end]]

        truth_bt = None
        inf_bt = None

        # --- build BedTools ---
        if not truth_s.empty:
            truth_bt = pybedtools.BedTool.from_dataframe(truth_s).sort().merge()
            t_true = sum(i.stop - i.start for i in truth_bt)
        else:
            t_true = 0

        if not inf_s.empty:
            inf_bt = pybedtools.BedTool.from_dataframe(inf_s).sort().merge()
            t_inf = sum(i.stop - i.start for i in inf_bt)
        else:
            t_inf = 0

        # --- clip if needed ---
        if sequence_length is not None:
            chroms = set()

            if truth_bt:
                chroms.update(i.chrom for i in truth_bt)
            if inf_bt:
                chroms.update(i.chrom for i in inf_bt)

            if chroms:
                clip = pybedtools.BedTool(
                    "\n".join(f"{c}\t0\t{sequence_length}" for c in chroms),
                    from_string=True
                )

                if truth_bt:
                    truth_bt = truth_bt.intersect(clip).sort().merge()
                    t_true = sum(i.stop - i.start for i in truth_bt)

                if inf_bt:
                    inf_bt = inf_bt.intersect(clip).sort().merge()
                    t_inf = sum(i.stop - i.start for i in inf_bt)

        # --- TP ---
        if truth_bt and inf_bt:
            tp = sum(i.stop - i.start for i in inf_bt.intersect(truth_bt).sort().merge())
        else:
            tp = 0

        # --- FP ---
        if inf_bt and truth_bt:
            fp = sum(i.stop - i.start for i in inf_bt.subtract(truth_bt).sort().merge())
        else:
            fp = t_inf

        # --- FN ---
        if truth_bt and inf_bt:
            fn = sum(i.stop - i.start for i in truth_bt.subtract(inf_bt).sort().merge())
        else:
            fn = t_true

        total_tp += tp
        total_fp += fp
        total_fn += fn
        total_true += t_true
        total_inf += t_inf

        if sequence_length is not None:
            total_bases += sequence_length

    result = [total_tp, total_fp, total_fn]

    # --- optional TN ---
    if return_tn and sequence_length is not None:
        tn = total_bases - (total_tp + total_fp + total_fn)
        result.append(tn)

    # --- optional precision/recall ---
    if return_precision_recall:
        precision = (total_tp / total_inf * 100) if total_inf else np.nan
        recall = (total_tp / total_true * 100) if total_true else np.nan
        result.extend([precision, recall])

    return tuple(result)