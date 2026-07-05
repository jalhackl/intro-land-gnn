import os
import subprocess
import sys
import shutil
import pandas as pd
import glob
import re
from resolve_indices_relate import build_trace_indices

def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


def load_idx(idx_file):
    with open(idx_file) as f:
        return [x.strip() for x in f.read().split(",") if x.strip()]


def _run_pipeline(
    trees,
    ref_list,
    tgt_list,
    idx_file,
    out_base,
    base_dir,
    prefix="trace",
    t=20000,
    chrom="1",
    posterior_threshold=0.9,
    physical_length_threshold=50000,
    genetic_distance_threshold=0.05
):
    """
    Full trace workflow
    """

    # Ensure TRACE-compatible tree filename
    if trees.endswith(".ts"):
        trees_copy = trees[:-3] + ".trees"

        if not os.path.exists(trees_copy):
            shutil.copy2(trees, trees_copy)

        trees = trees_copy

    npz_file = out_base + ".npz"

    summary_full_prefix = out_base + ".summary.full"
    summary_dir = out_base + ".summary.individual"

    os.makedirs(out_base, exist_ok=True)
    os.makedirs(base_dir, exist_ok=True)
    os.makedirs(summary_dir, exist_ok=True)

    # ------------------------------------------------
    # 1. build indices
    # ------------------------------------------------
    _, tgt_idx, tgt_id_map, tgt_name_to_index = build_trace_indices(
        ref_list, tgt_list, ploidy=2
    )

    with open(idx_file, "w") as f:
        f.write(",".join(map(str, tgt_idx)))

    inds = [str(i) for i in tgt_idx]

    # ------------------------------------------------
    # 2. trace-extract
    # ------------------------------------------------
    run(
        f"trace-extract "
        f"-f {trees} "
        f"-t {t} "
        f"-i $(cat {idx_file}) "
        f"-o {out_base}"
    )

    # ------------------------------------------------
    # 3. trace-infer
    # ------------------------------------------------
    infer_outputs = []

    for ind in inds:

        infer_prefix = os.path.join(
            #base_dir,
            out_base,
            f"infer.{prefix}.ind{ind}.xss.npz"
        )

        run(
            f"trace-infer "
            f"-i {ind} "
            f"--npz-files {npz_file} "
            f"--chroms {chrom} "
            f"-o {infer_prefix}"
        )

        infer_file = f"{infer_prefix}.{chrom}.xss.npz"

        if not os.path.exists(infer_file):
            raise FileNotFoundError(
                f"Expected TRACE output not found:\n{infer_file}"
            )

        infer_outputs.append(infer_file)

    # ============================================================
    # 4 + 5. trace-summarize + merge
    # ============================================================

    # Osingle value
    if isinstance(posterior_threshold, (int, float)):

        for ind, inf_file in zip(inds, infer_outputs):

            run(
                f"trace-summarize "
                f"-f {inf_file} "
                f"-c {chrom} "
                f"--posterior-threshold {posterior_threshold} "
                f"--physical-length-threshold {physical_length_threshold} "
                f"--genetic-distance-threshold {genetic_distance_threshold} "
                f"-o {os.path.join(summary_dir, f'ind{ind}.summary')}"
            )

        merge_summaries(
            summary_dir=summary_dir,
            output_file=out_base + ".ALL_merged_summary.txt",
            id_to_name=tgt_id_map
        )

    # list of thresholds

    else:

        for p_thr in posterior_threshold:

            p_thr_tag = f"threshold{p_thr}"

            # summarize per individual
            for ind, inf_file in zip(inds, infer_outputs):

                run(
                    f"trace-summarize "
                    f"-f {inf_file} "
                    f"-c {chrom} "
                    f"--posterior-threshold {p_thr} "
                    f"--physical-length-threshold {physical_length_threshold} "
                    f"--genetic-distance-threshold {genetic_distance_threshold} "
                    f"-o {os.path.join(summary_dir, f'ind{ind}.summary.{p_thr_tag}')}"
                )

            # merge per threshold
            merge_summaries(
                summary_dir=summary_dir,
                output_file=out_base + f".{p_thr_tag}.ALL_merged_summary.txt",
                id_to_name=tgt_id_map,
                file_middle_part= f"{p_thr_tag}.summary.txt",
                #file_suffix=p_thr_tag,
                export_bed=True,
                bed_output_file=out_base + f".{p_thr_tag}.inferred.tracts.bed"
            )


    done_file = out_base + ".done"
    with open(done_file, "w") as f:
        f.write("OK\n")


def main(
    trees,
    ref_list,
    tgt_list,
    scenario,
    rep,
    prefix,
    out_base,
    curr_folder,
    t=20000,
    log_file=None,
    chrom="1",
    posterior_threshold=0.9,
    physical_length_threshold=50000,
    genetic_distance_threshold=0.05
):

    base_dir = os.path.join(curr_folder, scenario, rep)

    idx_file = os.path.join(
        #base_dir,
        out_base,
        f"{prefix}.{rep}.tracehmm.idx"
    )

    try:
        _run_pipeline(
            trees=trees,
            ref_list=ref_list,
            tgt_list=tgt_list,
            idx_file=idx_file,
            out_base=out_base,
            base_dir=base_dir,
            prefix=prefix,
            t=t,
            chrom=chrom,
            posterior_threshold=posterior_threshold,
            physical_length_threshold=physical_length_threshold,
            genetic_distance_threshold=genetic_distance_threshold,
        )

    except Exception as e:
        print(f"TRACE pipeline failed: {e}")

        write_empty_outputs(out_base, posterior_threshold)

        # optionally write a .failed marker
        with open(out_base + ".failed", "w") as f:
            f.write(str(e) + "\n")


def write_empty_outputs(out_base, posterior_threshold):

    thresholds = (
        [posterior_threshold]
        if isinstance(posterior_threshold, (int, float))
        else posterior_threshold
    )

    for thr in thresholds:

        if isinstance(posterior_threshold, (int, float)):
            bed = out_base + ".inferred.tracts.bed"
            summary = out_base + ".ALL_merged_summary.txt"
        else:
            tag = f"threshold{thr}"
            bed = out_base + f".{tag}.inferred.tracts.bed"
            summary = out_base + f".{tag}.ALL_merged_summary.txt"

        # empty BED with header
        with open(bed, "w") as f:
            #f.write("chrom\tstart\tend\tindividual\n")
            f.write("")

        # empty summary
        with open(summary, "w") as f:
            f.write("")


def run_simple(
    trees,
    ref_list,
    tgt_list,
    out_folder,
    t=20000,
    chrom="1",
    posterior_threshold=0.9,
    physical_length_threshold=50000,
    genetic_distance_threshold=0.05
):
    """
    Call trace workflow (without specific names for the folders, only indicating, trees, ind lists, and out_folder)
    ----------
    trees : str
        Tree sequence file.

    ref_list : str
        Reference samples.

    tgt_list : str
        Target samples.

    out_folder : str
        Output directory.

    t : int
        Time parameter for trace.

    chrom : str
        Chromosome name.
    """

    os.makedirs(out_folder, exist_ok=True)

    idx_file = os.path.join(out_folder, "trace.idx")
    out_base = os.path.join(out_folder, "trace")

    _run_pipeline(
        trees=trees,
        ref_list=ref_list,
        tgt_list=tgt_list,
        idx_file=idx_file,
        out_base=out_base,
        base_dir=out_folder,
        t=t,
        chrom=chrom,
        posterior_threshold=posterior_threshold,
        physical_length_threshold=physical_length_threshold,
        genetic_distance_threshold=genetic_distance_threshold
    )



def merge_summaries(
    summary_dir,
    output_file,
    file_middle_part="summary.txt",
    file_suffix="",
    id_to_name=None,
    export_bed=True,
    bed_output_file=None
):

    # ------------------------------------------------
    # select files (old vs new naming)
    # ------------------------------------------------
    if file_suffix:
        pattern = os.path.join(summary_dir, f"ind*.{file_middle_part}.{file_suffix}")
    else:
        pattern = os.path.join(summary_dir, f"ind*.{file_middle_part}")

    files = sorted(glob.glob(pattern))

    if len(files) == 0:
        raise FileNotFoundError(f"No summary files found in {summary_dir}")

    all_dfs = []

    for path in files:
        base = os.path.basename(path)

        match = re.search(r"ind(\d+)", base)
        if not match:
            raise ValueError(f"Cannot parse individual id from {base}")

        individual_id = int(match.group(1))

        df = pd.read_csv(path, sep="\t")

        df["individual_id"] = individual_id

        if id_to_name is not None:
            df["individual"] = id_to_name[str(individual_id)]

        all_dfs.append(df)

    merged = pd.concat(all_dfs, ignore_index=True)
    merged.to_csv(output_file, sep="\t", index=False)

    # ------------------------------------------------
    # BED export
    # ------------------------------------------------
    if export_bed:
        if bed_output_file is None:
            bed_output_file = output_file.replace(
                "ALL_merged_summary.txt",
                "inferred.tracts.bed"
            )

        bed_df = merged[["chromosome", "start", "end", "individual"]].copy()
        bed_df.columns = ["chrom", "start", "end", "individual"]

        bed_df.to_csv(bed_output_file, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main(*sys.argv[1:])