import os
import glob
scenarios = main_config["scenarios"]
thresholds = trace_config["thresholds"]  # kept only for Python, NOT DAG splitting


# -------------------------
# TRACEHMM PIPELINE
# -------------------------
rule tracehmm_pipeline:
    input:
        trees=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.relate.trees"
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
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tracehmm.done"
        )

    params:
        out_base=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            wc.rep,
            f"{wc.prefix}.{wc.rep}.tracehmm.extract"
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
            "{prefix}.{rep}.tracehmm.log"
        )

    conda:
        "../../envs/trace-env.yaml"

    script:
        "run_trace_pipeline_smk.py"

# -------------------------
# ACCURACY PER THRESHOLD
# (threshold is ONLY analysis dimension now)
# -------------------------


import re

def extract_rep(path):
    m = re.search(r"\.(\d+)\.tracehmm\.extract|(\d+)\.true\.tracts", path)
    # fallback: extract any number near rep position
    m = re.search(r"\.(\d+)\.", path)
    return int(m.group(1)) if m else -1


rule tracehmm_accuracy:
    input:
        #rules.tracehmm_pipeline.output,
        true_tracts=lambda wc: [
            os.path.join(
                main_config["test_data_folder"],
                wc.scenario,
                str(rep),
                f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.true.tracts.bed"
            )
            for rep in range(scenario_configs[wc.scenario]["nrep"])
        ],

        inferred_tracts=lambda wc: glob.glob(
            os.path.join(
                main_config["test_data_folder"],
                wc.scenario,
                "*",
                f"{scenario_configs[wc.scenario]['output_prefix']}.*"
                f".tracehmm.extract.threshold{wc.threshold}.inferred.tracts.bed"
            )
        )

    output:
        accuracy=os.path.join(
            results_folder,
            "tracehmm",
            "{scenario}",
            "{prefix}.tracehmm.{scenario}.{threshold}.accuracy"
        ),

        metrics=os.path.join(results_folder, "tracehmm", "{scenario}", "{prefix}.tracehmm.{scenario}.{threshold}.metrics")




    run:

        print("all input pairs")
        print(input.true_tracts)
        print(input.inferred_tracts)

        true_map = {extract_rep(p): p for p in input.true_tracts}
        inf_map  = {extract_rep(p): p for p in input.inferred_tracts}

        common_reps = sorted(set(true_map) & set(inf_map))

        true_sorted = [true_map[r] for r in common_reps]
        inf_sorted  = [inf_map[r] for r in common_reps]


        precision, recall, tp, fp, fn, tn = cal_metrics_multiple(
            #input.true_tracts,
            #input.inferred_tracts,
            true_sorted,
            inf_sorted,

            scenario_configs[wildcards.scenario]["seq_len"]
        )

        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

        with open(output.accuracy, "w") as f:
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.threshold}\t"
                f"{precision}\t"
                f"{recall}\t"
                f"{f1}\n"
            )

        
                # Write detailed metrics including F1-score
        with open(output.metrics, "w") as f:
            f.write(
               # "scenario\tthreshold\tTP\tFP\tFN\tTN\tprecision\trecall\tF1\n"
               f"{wildcards.scenario}\t"
               f"{wildcards.threshold}\t"
               f"{tp}\t"
               f"{fp}\t"
               f"{fn}\t"
               f"{tn}\t"
               f"{precision}\t"
               f"{recall}\t"
               f"{f1}\n"
           )


# -------------------------
# SUMMARY
# -------------------------
rule tracehmm_accuracy_summary:
    input:
        expand(
            os.path.join(
                results_folder,
                "tracehmm",
                "{scenario}",
                "{prefix}.tracehmm.{scenario}.{threshold}.accuracy"
            ),
            scenario=scenarios,
            threshold=thresholds,
            prefix=[scenario_configs[s]["output_prefix"] for s in scenario_configs],
        )

    output:
        os.path.join(
            results_folder,
            "tracehmm",
            "tracehmm_accuracy_aggregated.txt"
        )

    shell:
        """
        cat {input} | sed '1i scenario\tthreshold\tprecision\trecall\tf1' > {output}
        """


#single tract file processing
rule tracehmm_rep_threshold_metrics:
    input:
        rules.tracehmm_pipeline.output,
        true=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.rep}.true.tracts.bed"
        ),



    output:
        metrics=os.path.join(
            results_folder,
            "tracehmm",
            "{scenario}",
            "{prefix}.rep_{rep}",
            "threshold_{threshold}",
            "metrics_ind.txt"
        )

    params:
        seq_len=lambda wc: scenario_configs[wc.scenario]["seq_len"],
        inferred=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.rep}.tracehmm.extract.threshold{wc.threshold}.inferred.tracts.bed"
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
            f.write(
                "scenario\trep\tthreshold\tTP\tFP\tFN\tTN\tprecision\trecall\n"
            )
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.rep}\t"
                f"{wildcards.threshold}\t"
                f"{tp}\t{fp}\t{fn}\t{tn}\t"
                f"{precision}\t{recall}\n"
            )



rule tracehmm_metrics_collect:
    input:
        [
            os.path.join(
                results_folder,
                "tracehmm",
                scenario,
                f"{scenario_configs[scenario]['output_prefix']}.rep_{rep}",
                f"threshold_{threshold}",
                "metrics_ind.txt"
            )
            for scenario in scenarios
            for rep in range(scenario_configs[scenario]["nrep"])
            for threshold in thresholds
        ]
    output:
        os.path.join(
            results_folder,
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



rule tracehmm_metrics_summary:
    input:
        os.path.join(
            results_folder,
            "tracehmm",
            "tracehmm_metrics_all.csv"
        )

    output:
        os.path.join(
            results_folder,
            "tracehmm",
            "tracehmm_metrics_summary.csv"
        )

    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")

        summary = (
            df.groupby(["scenario", "threshold"], as_index=False)
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

