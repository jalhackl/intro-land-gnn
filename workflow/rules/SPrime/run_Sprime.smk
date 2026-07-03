import os
from sprime_additional_functions import create_map_file, process_sprime_output

from sprime_mapping import map_sprime_segments

sys.path.append("workflow/scripts")

from evaluate_utils import cal_accuracy, cal_accuracy_multiple, cal_metrics_multiple

data_folder = main_config["test_data_folder"]
results_folder = main_config["results_folder"]
scenario_list = list(scenario_configs.keys())
threshold_list = sprime_config["threshold_list"]

# get replicates for a scenario after checkpoint finishes
def get_replicates_for_scenario(scenario):
    ck = checkpoints.simulate_gaishisim.get(scenario=scenario)
    scenario_dir = os.path.join(data_folder, scenario)
    replicates = sorted([
        d for d in os.listdir(scenario_dir)
        if os.path.isdir(os.path.join(scenario_dir, d)) and not d.startswith(".")
    ])
    return replicates

# Rule: Run Sprime on each scenario/replicate/threshold
rule sprime_run:
    input:
        vcf=lambda wc: checkpoints.get_biallelic_vcf.get(
            scenario=wc.scenario,
            replicate=wc.replicate,
            output_prefix=f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}"
        ).output.vcf,
        ref_list=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            str(wc.replicate),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.ref.ind.list"
        )
    output:
        score=os.path.join(
    results_folder, "sprime", "{scenario}", "{replicate}", "sprime.1src.out.{threshold}.score"
    )
    params:
        sprime_exec=sprime_config["sprime_exec"],
        threshold=lambda wc: float(wc.threshold),
        output_prefix=os.path.join(results_folder, "sprime", "{scenario}", "{replicate}", "sprime.1src.out.{threshold}"),
        mut_rate=lambda wc: float(scenario_configs[wc.scenario]['mut_rate']),
        rec_rate=lambda wc: float(scenario_configs[wc.scenario]['rec_rate']),
        seq_len=lambda wc: int(scenario_configs[wc.scenario]['seq_len']),
    threads: 1
    resources:
        time_min=60,
        mem_mb=2000,
        cpus=1
    benchmark:
        "benchmarks/sprime/{scenario}/{replicate}/{threshold}.txt"
    run:
        output_dir = os.path.dirname(params.output_prefix)
        os.makedirs(output_dir, exist_ok=True)

        # Unique map file per run (probably not necessary, only per scenario)
        map_file_config = sprime_config.get("map_file", None)
        if map_file_config is None:
            recomb_cM = params.seq_len * params.rec_rate * 100 
            map_file = os.path.join(output_dir, f"sim_{wildcards.scenario}_{wildcards.replicate}.map")
            create_map_file(recomb_cM, params.seq_len, map_file)
        else:
            map_file = map_file_config

        shell(
            "java -Xmx2g -jar {params.sprime_exec} "
            "gt={input.vcf} "
            "outgroup={input.ref_list} "
            "map={map_file:q} "
            "out={params.output_prefix} "
            "minscore={params.threshold} "
            "mu={params.mut_rate}"
        )


rule sprime_map_results:
    input:
        vcf=lambda wc: checkpoints.get_biallelic_vcf.get(
            scenario=wc.scenario,
            replicate=wc.replicate,
            output_prefix=f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}"
        ).output.vcf,
        tgt_list=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            str(wc.replicate),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.tgt.ind.list"
        ),
        scores=rules.sprime_run.output.score
    output:
        inferred_tracts=os.path.join(results_folder, "sprime", "{scenario}", "{replicate}", "sprime.1src.out.{threshold}.bed"),
    params:
        merge_distance=0
    threads: 1
    resources:
        time_min=60,
        mem_mb=2000,
        cpus=1
    benchmark:
        "benchmarks/sprime/{scenario}/{replicate}_{threshold}_maptracts.txt"
    run:
        with open(input.tgt_list, "r") as f:
            tgt_inds = [line.strip() for line in f]

        map_sprime_segments(
        input.scores,
        input.vcf,
        out_file=output.inferred_tracts,
        target_individuals=None,
        segment_fraction=None,
        min_snps=None,
        phased=True,
        merge_distance=params.merge_distance,
        only_tract_output=True,
        return_full_records=False,
        single_snp_bed=False, 
        )
                


#Rule: Process mapped Sprime output to generate accuracy
rule sprime_process_mapped_output:
    input:
        inferred_tracts=rules.sprime_map_results.output.inferred_tracts,
        true_tracts=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.true.tracts.bed"
        )
    output:
        accuracy=os.path.join(results_folder, "sprime", "{scenario}", "{replicate}", "sprime.1src.out.{threshold}.accuracy")
    threads: 1
    resources:
        time_min=60,
        mem_mb=2000,
        cpus=1
    run:
        precision, recall = cal_accuracy_samples(input.true_tracts, input.inferred_tracts)
        with open(output.accuracy, 'w') as f:
            f.write(f"{wildcards.scenario}\t{wildcards.replicate}\t{wildcards.threshold}\t{precision}\t{recall}\n")

#Function to get all accuracy files 
def get_all_accuracy_files(wildcards):
    all_files = []
    for scenario in scenario_list:
        replicates = get_replicates_for_scenario(scenario)
        for rep in replicates:
            for threshold in threshold_list:
                all_files.append(os.path.join(
                 results_folder, "sprime", scenario, rep, f"sprime.1src.out.{threshold}.accuracy"
                ))
    return all_files

#Rule: Aggregate all accuracy results into one summary table
rule sprime_accuracy:
    input:
        accuracy_files=get_all_accuracy_files
    output:
        accuracy_table=os.path.join(results_folder, "sprime", "sprime_accuracy.txt")
    shell:
        """
        cat {input.accuracy_files} | sed '1i scenario\treplicate\tcutoff\tprecision\trecall' > {output.accuracy_table}
        """

#trigger all sprime runs & processing
rule sprime_all:
    input:
        get_all_accuracy_files



rule sprime_accuracy_aggregated:
    input:
        true_tracts=lambda wc: [
            os.path.join(
                data_folder,
                wc.scenario,
                str(rep),
                f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.true.tracts.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ],
        inferred_tracts=lambda wc: [
            os.path.join(
                results_folder, "sprime", wc.scenario, rep,
                f"sprime.1src.out.{wc.threshold}.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ]
    params:
        seq_len=lambda wc: int(scenario_configs[wc.scenario]["seq_len"])
    output:
        out=os.path.join(results_folder, "sprime", "{scenario}", "sprime.aggregated.{threshold}.txt"),
        metrics=os.path.join(results_folder, "sprime", "{scenario}", "sprime.aggregated.{threshold}.metrics")
    run:
        print("Scenario:", wildcards.scenario)
        print("Threshold:", wildcards.threshold)
        print("seq_len:", params.seq_len)
        print("True tracts:", input.true_tracts)
        print("Inferred tracts:", input.inferred_tracts)


        (
            precision,
            recall,
            total_true_positive,
            total_false_positive,
            total_false_negative,
            total_true_negative
        ) = cal_metrics_multiple(
            input.true_tracts,
            input.inferred_tracts,
            params.seq_len
        )

        # Compute F1-score safely (avoid division by zero)
        if precision + recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0

        # Write precision/recall summary
        with open(output.out, "w") as f:
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.threshold}\t"
                f"{precision}\t"
                f"{recall}\n"
            )

        # Write detailed metrics including F1-score
        with open(output.metrics, "w") as f:
            f.write(
               # "scenario\tthreshold\tTP\tFP\tFN\tTN\tprecision\trecall\tF1\n"
               f"{wildcards.scenario}\t"
               f"{wildcards.threshold}\t"
               f"{total_true_positive}\t"
               f"{total_false_positive}\t"
               f"{total_false_negative}\t"
               f"{total_true_negative}\t"
               f"{precision}\t"
               f"{recall}\t"
               f"{f1}\n"
           )




rule sprime_accuracy_summary:
    input:
        expand(os.path.join(results_folder, "sprime", "{scenario}", "sprime.aggregated.{threshold}.txt"),
               scenario=scenario_list, threshold=threshold_list)
    output:
        os.path.join(results_folder, "sprime", "sprime_accuracy_aggregated.txt")
    shell:
        """
        cat {input} | sed '1i scenario\tthreshold\tprecision\trecall' \
            > {output}
        """



rule sprime_metrics_summary:
    input:
        expand(os.path.join(results_folder, "sprime", "{scenario}", "sprime.aggregated.{threshold}.metrics"),
               scenario=scenario_list, threshold=threshold_list)
    output:
        os.path.join(results_folder, "sprime", "sprime_metrics_aggregated.txt")
    shell:
        """
        cat {input} | sed '1i scenario\tthreshold\tTP\tFP\tFN\tTN\tF1' \
            > {output}
        """


