import os
import glob
import re


sys.path.append("workflow/scripts")
from evaluate_utils import cal_accuracy, cal_accuracy_multiple, cal_metrics_multiple, cal_metrics_single


scenarios = main_config["scenarios"]
thresholds = trace_config["thresholds"]

tree_types = [
    ".relate.trees",
    ".ts",
    "_tsinfer.ts"
]

# -------------------------
# TRACEHMM PIPELINE
# -------------------------
rule tracehmm_pipeline:
    input:
        trees=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}{tree_type}"
        ),
        ref_list=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.ref.ind.list"
        ),
        tgt_list=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tgt.ind.list"
        )

    output:
        done=os.path.join(
            main_config["results_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.{tree_type}.tracehmm.done"
        )

    params:
        out_base=lambda wc: os.path.join(
            main_config["results_folder"],
            wc.scenario,
            wc.rep,
            f"{wc.prefix}.{wc.rep}.{wc.tree_type}.tracehmm.extract"
        ),
        curr_folder=main_config["test_data_folder"],
        thresholds=trace_config["thresholds"],
        t=trace_config["t"],
        chrom="1",
        physical_length_threshold=trace_config["physical_length_threshold"],
        genetic_distance_threshold=trace_config["genetic_distance_threshold"]

    log:
        os.path.join(
            "logs",
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.{tree_type}.tracehmm.log"
        )

    conda:
        "../../envs/trace-env.yaml"

    script:
        "run_trace_pipeline_smk.py"



# -------------------------
#  METRICS
# -------------------------
rule tracehmm_rep_threshold_metrics:
    input:
        rules.tracehmm_pipeline.output,
        true=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.rep}.true.tracts.bed"
        )

    output:
        metrics=os.path.join(
            main_config["results_folder"],
            "tracehmm",
            "{scenario}",
            "{tree_type}",
            "{prefix}.rep_{rep}",
            "threshold_{threshold}",
            "metrics_ind.txt"
        )

    params:
        seq_len=lambda wc: scenario_configs[wc.scenario]["seq_len"],
        inferred=lambda wc: os.path.join(
            main_config["results_folder"],
            wc.scenario,
            str(wc.rep),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.rep}.{wc.tree_type}.tracehmm.extract.threshold{wc.threshold}.inferred.tracts.bed"
        )

    run:

        tp, fp, fn, tn = cal_metrics_single(
            input.true,
            params.inferred,
            sequence_length=params.seq_len,
            return_tn=True
        )

        precision = tp / (tp + fp) if (tp + fp) else 0
        recall = tp / (tp + fn) if (tp + fn) else 0

        with open(output.metrics, "w") as f:
            '''
            f.write(
                "scenario\trep\tthreshold\ttree_type\tTP\tFP\tFN\tTN\tprecision\trecall\n"
            )
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.rep}\t"
                f"{wildcards.threshold}\t"
                f"{wildcards.tree_type}\t"
                f"{tp}\t{fp}\t{fn}\t{tn}\t"
                f"{precision}\t{recall}\n"
            )
            '''
            #including tree_type
            f.write(
                "scenario\ttree_type\trep\tthreshold\tTP\tFP\tFN\tTN\tprecision\trecall\n"
            )
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.tree_type}\t"
                f"{wildcards.rep}\t"
                f"{wildcards.threshold}\t"
                f"{tp}\t{fp}\t{fn}\t{tn}\t"
                f"{precision}\t{recall}\n"
            )

# COLLECT METRICS
rule tracehmm_metrics_collect:
    input:
        [
            os.path.join(
                main_config["results_folder"],
                "tracehmm",
                scenario,
                tree_type,
                f"{scenario_configs[scenario]['output_prefix']}.rep_{rep}",
                f"threshold_{threshold}",
                "metrics_ind.txt"
            )
            for scenario in scenarios
            for tree_type in tree_types
            for rep in range(scenario_configs[scenario]["nrep"])
            for threshold in thresholds
        ]

    output:
        os.path.join(
            main_config["results_folder"],
            "tracehmm",
            "tracehmm_metrics_all.csv"
        )

    run:
        import pandas as pd

        dfs = [pd.read_csv(f, sep="\t") for f in input]

        pd.concat(dfs, ignore_index=True).to_csv(
            output[0],
            sep="\t",
            index=False
        )

# SUMMARY
rule tracehmm_metrics_summary:
    input:
        os.path.join(
            main_config["results_folder"],
            "tracehmm",
            "tracehmm_metrics_all.csv"
        )

    output:
        os.path.join(
            main_config["results_folder"],
            "tracehmm",
            "tracehmm_metrics_summary.csv"
        )

    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")

        summary = (
            df.groupby(["scenario", "tree_type", "threshold"], as_index=False)
              .agg({
                  "TP": "sum",
                  "FP": "sum",
                  "FN": "sum",
                  "TN": "sum"
              })
        )

        summary["precision"] = (
            summary["TP"] /
            (summary["TP"] + summary["FP"]) * 100
        )

        summary["recall"] = (
            summary["TP"] /
            (summary["TP"] + summary["FN"]) * 100
        )

        summary["F1"] = (
            2 * summary["precision"] * summary["recall"] /
            (summary["precision"] + summary["recall"])
        )

        summary.to_csv(
            output[0],
            sep="\t",
            index=False
        )