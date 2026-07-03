import numpy as np
import os

from hmmix_additional_functions_v082 import *

import make_mutationrate
from make_mutationrate import make_mutation_rate
from hmm_functions import HMMParam, write_HMM_to_file

sys.path.append("workflow/scripts")
from evaluate_utils import *


# Configuration from SPrime workflow
data_folder = main_config["test_data_folder"]
results_folder = main_config["results_folder"]
scenario_list = list(scenario_configs.keys())
cutoff_list = hmmix_config.get("cutoff_list", [0.5])  # analogous to thresholds

try:
    binary = hmmix_config["binary"]
except KeyError:
    binary = False

if binary:
    ref_set = ["0", "1"]
else:
    ref_set = ["A", "C", "G", "T"]

ruleorder: skov_run > hmmix_accuracy_aggregated

# Get replicates (same as SPrime)
def get_replicates_for_scenario(scenario):
    scenario_dir = os.path.join(data_folder, scenario)
    replicates = sorted([d for d in os.listdir(scenario_dir)
                         if os.path.isdir(os.path.join(scenario_dir, d)) and not d.startswith(".")])
    return replicates


rule skov_run:
    input:
        vcf=lambda wc: checkpoints.get_biallelic_vcf.get(
            scenario=wc.scenario,
            replicate=wc.replicate,
            output_prefix=f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}"
        ).output.vcf,
        ref_list=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.ref.ind.list"
        ),
        tgt_list=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.tgt.ind.list"
        ),
        true_tracts=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.true.tracts.bed"
        )
        
    output:
        prob_file=os.path.join(results_folder, "hmmix_classic", "{scenario}", "{replicate}", "probabilities.txt"),
   
    
    params:
        out_folder=lambda wc: os.path.join(results_folder, "hmmix_classic", wc.scenario, wc.replicate),
        output_prefix=lambda wc: scenario_configs[wc.scenario]['output_prefix'],
        src_id=lambda wc: scenario_configs[wc.scenario]['src_id'],
        ref_id=lambda wc: scenario_configs[wc.scenario]['ref_id']
    threads: 16
    resources:
        mem_gb=32
    run:
        os.makedirs(params.out_folder, exist_ok=True)

        # Load reference and target individuals
        with open(input.ref_list, 'r') as f:
            ref_lines = [line.strip() for line in f]
        with open(input.tgt_list, 'r') as f:
            tgt_lines = [line.strip() for line in f]

        # Outgroup and ingroup processing
        make_out_group_custom_mut(ref_lines, None, [input.vcf], os.path.join(params.out_folder, "output_ref_vs.txt"), [None], [None], ref_set=ref_set)

        make_ingroup_custom_mut(tgt_lines, None, [input.vcf], os.path.join(params.out_folder, "output_tgt_vs"), os.path.join(params.out_folder, "output_ref_vs.txt"), [None], ref_set=ref_set)

        # Mutation rate
        make_mutation_rate(os.path.join(params.out_folder, "output_ref_vs.txt"), os.path.join(params.out_folder, "mutrates_outgroup1.out"), None, 100000)

        # Initial HMM guess
        new_HMM = HMMParam([params.ref_id, params.src_id], [0.5, 0.5], [[0.99, 0.01], [0.02, 0.98]], [0.03, 0.3])
        write_HMM_to_file(new_HMM, os.path.join(params.out_folder, "hmm_guesses.json"))

        # Train and decode
        infiles, trained_files = train_hmm_individuals(tgt_lines, os.path.join(params.out_folder, "hmm_guesses.json"),
                                                        os.path.join(params.out_folder, "mutrates_outgroup1.out"),
                                                        out_folder=params.out_folder, window_size=1000, haploid=False)

        
        #version 0.82 needs also chromsome
        all_segments = decode_hmm_individuals(tgt_lines, trained_files, os.path.join(params.out_folder, "mutrates_outgroup1.out"),1,
                                              out_folder=params.out_folder, window_size=1000, haploid=False)

        # old version
        #all_segments = decode_hmm_individuals(tgt_lines, trained_files, os.path.join(params.out_folder, "mutrates_outgroup1.out"),
        #                                      out_folder=params.out_folder, window_size=1000, haploid=False)
        all_segments.to_csv(output.prob_file, index=False)

        # Process output for each cutoff
        inferred_tract_files = process_output(all_segments, params.out_folder, params.output_prefix, params.src_id, cutoff_list=cutoff_list, return_filenames=True)
        for inferred_tracts in inferred_tract_files:
            cutoff = inferred_tracts[0]
            precision, recall = cal_accuracy_samples(input.true_tracts, inferred_tracts[1])
            acc_file = inferred_tracts[1].rsplit('.', 1)[0] + ".accuracy"
            with open(acc_file, 'w') as f:
                f.write(f"{wildcards.scenario}\t{wildcards.replicate}\t{cutoff}\t{precision}\t{recall}\n")

# Function to get all probability files for a given scenario
def get_skov_files_for_scenario(wc):
    replicates = get_replicates_for_scenario(wc.scenario)
    return expand(
        os.path.join(
            results_folder,
            "hmmix_classic",
            wc.scenario,
            "{replicate}",
            "probabilities.txt"
        ),
        replicate=replicates
    )

rule hmmix_scenario_done:
    input:
        get_skov_files_for_scenario
    output:
        os.path.join(
            results_folder,
            "hmmix_classic",
            "{scenario}",
            ".hmmix_done"
        )
    shell:
        "touch {output}"


rule hmmix_accuracy_aggregated:
    input:
        true_tracts=lambda wc: [
            os.path.join(
                data_folder,
                wc.scenario,
                str(rep),
                f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.true.tracts.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ]

    params:
        seq_len=lambda wc: int(scenario_configs[wc.scenario]["seq_len"]),
        inferred_tracts=lambda wc: [
            os.path.join(
                results_folder,
                "hmmix_classic",
                wc.scenario,
                rep,
                #f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.src{scenario_configs[wc.scenario]['src_id']}.cutoff{wc.cutoff}.bed"
                f"{scenario_configs[wc.scenario]['output_prefix']}{wc.cutoff}.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ]
    output:
        out=os.path.join(
            results_folder,
            "hmmix_classic",
            "{scenario}",
            "hmmix_classic.aggregated.{cutoff}.txt"
        ),
        metrics=os.path.join(
            results_folder,
            "hmmix_classic",
            "{scenario}",
            "hmmix_classic.aggregated.{cutoff}.metrics"
        )
    run:
        (
            precision,
            recall,
            total_true_positive,
            total_false_positive,
            total_false_negative,
            total_true_negative
        ) = cal_metrics_multiple(
            input.true_tracts,
            params.inferred_tracts,
            params.seq_len
        )

        if precision + recall > 0:
            f1 = 2 * (precision * recall) / (precision + recall)
        else:
            f1 = 0.0

        with open(output.out, "w") as f:
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.cutoff}\t"
                f"{precision}\t"
                f"{recall}\n"
            )

        with open(output.metrics, "w") as f:
            f.write(
                f"{wildcards.scenario}\t"
                f"{wildcards.cutoff}\t"
                f"{total_true_positive}\t"
                f"{total_false_positive}\t"
                f"{total_false_negative}\t"
                f"{total_true_negative}\t"
                f"{precision}\t"
                f"{recall}\t"
                f"{f1}\n"
            )


rule hmmix_accuracy_summary:
    input:
        expand(
            os.path.join(
                results_folder,
                "hmmix_classic",
                "{scenario}",
                "hmmix_classic.aggregated.{cutoff}.txt"
            ),
            scenario=scenario_list,
            cutoff=cutoff_list
        )
    output:
        os.path.join(
            results_folder,
            "hmmix_classic",
            "hmmix_classic_accuracy_aggregated.txt"
        )
    shell:
        """
        cat {input} | sed '1i scenario\tcutoff\tprecision\trecall' \
            > {output}
        """


rule hmmix_metrics_summary:
    input:
        expand(
            os.path.join(
                results_folder,
                "hmmix_classic",
                "{scenario}",
                "hmmix_classic.aggregated.{cutoff}.metrics"
            ),
            scenario=scenario_list,
            cutoff=cutoff_list
        )
    output:
        os.path.join(
            results_folder,
            "hmmix_classic",
            "hmmix_classic_metrics_aggregated.txt"
        )
    shell:
        """
        cat {input} | sed '1i scenario\tcutoff\tTP\tFP\tFN\tTN\tprecision\trecall\tF1' \
            > {output}
        """