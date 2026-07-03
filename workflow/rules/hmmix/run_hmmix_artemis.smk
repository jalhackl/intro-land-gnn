import os
import numpy as np
from collections import defaultdict


from hmmix_additional_functions_v082 import *
from make_mutationrate import make_mutation_rate
from hmm_functions import HMMParam, write_HMM_to_file

# Configurations
data_folder = main_config["test_data_folder"]
results_folder = main_config["results_folder"]
scenario_list = list(scenario_configs.keys())
cutoff_list = hmmix_config.get("cutoff_list", [0.5])
binary = hmmix_config.get("binary", False)
ref_set = ["0", "1"] if binary else ["A", "C", "G", "T"]

skov_output_dir = os.path.join(results_folder, "hmmix_artemis")

default_alpha = hmmix_config.get("hybrid_default_value", 0)
artemis_window_size = hmmix_config.get("artemis_window_size", 1000)



#------------------------------------------------------------

def get_replicates_after_simulation(wc):
    ckpt = checkpoints.simulate_gaishisim.get(scenario=wc.scenario)
    scenario_dir = os.path.dirname(ckpt.output.done)

    replicates = sorted(
        d for d in os.listdir(scenario_dir)
        if os.path.isdir(os.path.join(scenario_dir, d)) and not d.startswith(".")
    )

    return replicates


def get_target_inds_after_simulation(wc, replicate):
    scenario_dir = os.path.join(data_folder, wc.scenario)

    tgt_file = os.path.join(
        scenario_dir,
        replicate,
        f"{scenario_configs[wc.scenario]['output_prefix']}.{replicate}.tgt.ind.list"
    )

    with open(tgt_file) as f:
        inds = [line.strip() for line in f]

    return inds


def get_all_hybrid_outputs(wc):
    ckpt = checkpoints.simulate_gaishisim.get(scenario=wc.scenario)
    scenario_dir = os.path.dirname(ckpt.output.done)

    outputs = []

    replicates = sorted(
        d for d in os.listdir(scenario_dir)
        if os.path.isdir(os.path.join(scenario_dir, d)) and not d.startswith(".")
    )

    for rep in replicates:
        tgt_file = os.path.join(
            scenario_dir,
            rep,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.tgt.ind.list"
        )

        with open(tgt_file) as f:
            inds = [line.strip() for line in f]

        for ind in inds:
            outputs.append(
                os.path.join(
                    skov_output_dir,
                    wc.scenario,
                    rep,
                    f"output_tgt_vs.{ind}.hybrid.decoded.txt.diploid.txt"
                )
            )

    return outputs



ruleorder:
    hmmix_artemis_run > hmmix_artemis_with_alpha >  hmmix_hybrid_decode > hmmix_all_done

#------------------------------------------------------------
# Rule: hmmix_artemis_run
rule hmmix_artemis_run:
    input:
        vcf=lambda wc: checkpoints.get_biallelic_vcf.get(
            scenario=wc.scenario,
            replicate=wc.replicate,
            output_prefix=f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}"
        ).output.vcf,
        ref_list=lambda wc: os.path.join(
            data_folder, wc.scenario, wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.ref.ind.list"
        ),
        tgt_list=lambda wc: os.path.join(
            data_folder, wc.scenario, wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.tgt.ind.list"
        ),
        true_tracts=lambda wc: os.path.join(
            data_folder, wc.scenario, wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.true.tracts.bed"
        )
    output:
        mutrates=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}", "mutrates_outgroup1.out"
        )
    params:
        out_folder=lambda wc: os.path.join(skov_output_dir, wc.scenario, wc.replicate),
        output_prefix=lambda wc: scenario_configs[wc.scenario]['output_prefix'],
        src_id=lambda wc: scenario_configs[wc.scenario]['src_id'],
        ref_id=lambda wc: scenario_configs[wc.scenario]['ref_id'],
        prob_file=lambda wc: os.path.join(skov_output_dir, wc.scenario, wc.replicate, "probabilities.txt")
    threads: 16
    resources:
        mem_gb=32
    run:
        os.makedirs(params.out_folder, exist_ok=True)

        with open(input.ref_list) as f:
            ref_lines = [line.strip() for line in f]
        with open(input.tgt_list) as f:
            tgt_lines = [line.strip() for line in f]

        make_out_group_custom_mut(
            ref_lines, None, [input.vcf],
            os.path.join(params.out_folder, "output_ref_vs.txt"),
            [None], [None], ref_set=ref_set
        )

        make_ingroup_custom_mut(
            tgt_lines, None, [input.vcf],
            os.path.join(params.out_folder, "output_tgt_vs"),
            os.path.join(params.out_folder, "output_ref_vs.txt"),
            [None], ref_set=ref_set
        )

        make_mutation_rate(
            os.path.join(params.out_folder, "output_ref_vs.txt"),
            output.mutrates, None, 100000
        )

        new_HMM = HMMParam(
            [params.ref_id, params.src_id],
            [0.5, 0.5],
            [[0.99, 0.01], [0.02, 0.98]],
            [0.03, 0.3]
        )

        write_HMM_to_file(new_HMM,
                          os.path.join(params.out_folder, "hmm_guesses.json"))

        infiles, trained_files = train_hmm_individuals(
            tgt_lines,
            os.path.join(params.out_folder, "hmm_guesses.json"),
            output.mutrates,
            out_folder=params.out_folder,
            window_size=1000,
            haploid=False
        )

        all_segments = decode_hmm_individuals(
            tgt_lines,
            trained_files,
            output.mutrates,
            1,
            out_folder=params.out_folder,
            window_size=1000,
            haploid=False
        )

        all_segments.to_csv(params.prob_file, index=False)

        inferred_tract_files = process_output(
            all_segments,
            params.out_folder,
            params.output_prefix,
            params.src_id,
            cutoff_list=cutoff_list,
            return_filenames=True
        )

        for cutoff, bed_file in inferred_tract_files:
            precision, recall = cal_accuracy(input.true_tracts, bed_file)
            acc_file = bed_file.rsplit('.', 1)[0] + ".accuracy"
            with open(acc_file, 'w') as f:
                f.write(f"{wildcards.scenario}\t{wildcards.replicate}\t{cutoff}\t{precision}\t{recall}\n")




rule hmmix_artemis_with_alpha:
    input:
        mutrates=lambda wc: os.path.join(
            skov_output_dir, wc.scenario, wc.replicate,
            "mutrates_outgroup1.out"
        )
    params:
        trained=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "output_tgt_vs.{ind}trained.json"
        ),
        obs=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "output_tgt_vs.{ind}.txt"
        ),
        art_out=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "artemis_tgt_vs.{ind}"
        ),
        art_plot=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "artemis_tgt_vs.{ind}.pdf"
        ),
    output:
        art_out=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "artemis_tgt_vs.{ind}"
        ),
        art_plot=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "artemis_tgt_vs.{ind}.pdf"
        ),
        alpha=os.path.join(
            skov_output_dir, "{scenario}", "{replicate}",
            "alpha_tgt_vs.{ind}"
        )
    conda:
        "../../envs/hmmix-env.yaml"
    shell:
        r"""
        export QT_QPA_PLATFORM='offscreen'
        set +e

        if [ ! -f "{params.trained}" ]; then
            echo "WARNING missing param {params.trained}, skipping"
            echo {default_alpha} > {output.alpha}
            touch {output.art_out} {output.art_plot}
            exit 0
        fi

        hmmix artemis \
            -param={params.trained} \
            -obs={params.obs} \
            -mutrates={input.mutrates} \
            -out={params.art_out} \
            -out_plot={params.art_plot} \
            -window_size={artemis_window_size}

        status=$?

        if [ $status -ne 0 ]; then
            echo "WARNING: Artemis failed. Writing default alpha."
            echo {default_alpha} > {output.alpha}
            touch {output.art_out} {output.art_plot}
        else
            Rscript --vanilla workflow/rules/hmmix/extract_artemis_alpha.R \
                {params.art_out} {output.alpha}
        fi

        exit 0
        """

#------------------------------------------------------------
# Rule: hybrid decode
rule hmmix_hybrid_decode:
    input:
        mut_rate=os.path.join(skov_output_dir, "{scenario}", "{replicate}", "mutrates_outgroup1.out"),
        alpha=os.path.join(skov_output_dir, "{scenario}", "{replicate}", "alpha_tgt_vs.{ind}")
    output:
        os.path.join(skov_output_dir, "{scenario}", "{replicate}", "output_tgt_vs.{ind}.hybrid.decoded.txt.diploid.txt")
    params:
        output_file = os.path.join(skov_output_dir, "{scenario}", "{replicate}", "output_tgt_vs.{ind}.hybrid.decoded.txt"),
        weights_bed=config.get("weights_bed", ""),
        obs=os.path.join(skov_output_dir, "{scenario}", "{replicate}", "output_tgt_vs.{ind}.txt"),
        trained=os.path.join(skov_output_dir, "{scenario}", "{replicate}", "output_tgt_vs.{ind}trained.json"),
    log:
        os.path.join(skov_output_dir, "{scenario}", "{replicate}", "logs", "hybriddecode.{ind}.log")
    resources:
        mem_gb=32,
        cpus=16
    conda:
        "../../envs/hmmix-env.yaml"
    threads: 1
    shell:
        """
        python workflow/rules/hmmix/hmmix_decode_wrapper.py \
            --obs {params.obs} \
            --mutrates {input.mut_rate} \
            --param {params.trained} \
            --out {params.output_file} \
            --alpha {input.alpha} \
            --weights "{params.weights_bed}"
        """

rule combine_hmmix_decoded:
    input:
        lambda wc: expand(
            os.path.join(
                skov_output_dir,
                wc.scenario,
                wc.replicate,
                "output_tgt_vs.{ind}.hybrid.decoded.txt.diploid.txt"
            ),
            ind=get_target_inds_after_simulation(wc, wc.replicate)
        )
    output:
        os.path.join(
            skov_output_dir,
            "{scenario}",
            "{replicate}",
            "output_tgt_vs.combined.decoded.txt.diploid.txt"
        )
    run:
        import pandas as pd
        import os
        import re

        dfs = []
        for file in input:
            sample = re.search(
                r'output_tgt_vs\.(.*?)\.hybrid',
                os.path.basename(file)
            ).group(1)

            df = pd.read_csv(file, sep="\t")

            if df.empty:
                print(f"Warning! Empty df in combine_hmmix_decoded: {file}")
                continue

            df["sample"] = sample
            dfs.append(df)

        pd.concat(dfs, ignore_index=True).to_csv(output[0], sep="\t", index=False)



rule hmmix_artemis_processing:
    input:
        decoded=os.path.join(
            skov_output_dir,
            "{scenario}",
            "{replicate}",
            "output_tgt_vs.combined.decoded.txt.diploid.txt"
        ),
        true_tracts=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            str(wc.replicate),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.true.tracts.bed"
        )
    output:
        prob_file=os.path.join(
            results_folder,
            "hmmix_artemis",
            "{scenario}",
            "{replicate}",
            "combined_segments.txt"
        )
    params:
        out_folder=lambda wc: os.path.join(results_folder, "hmmix_artemis", wc.scenario, str(wc.replicate)),
        output_prefix=lambda wc: scenario_configs[wc.scenario]['output_prefix'],
        src_id=lambda wc: scenario_configs[wc.scenario]['src_id']
    run:
        import os
        import shutil
        os.makedirs(params.out_folder, exist_ok=True)

        # Save combined decoded table (like probabilities)
        all_segments = pd.read_csv(input.decoded, sep="\t")

        # Process all cutoffs at once
        new_output_prefix = params.output_prefix + ".hybrid."
        inferred_tract_files = process_output(
            all_segments,
            params.out_folder,
            new_output_prefix,
            params.src_id,
            cutoff_list=cutoff_list,
            return_filenames=True
        )

        # Compute precision/recall for each cutoff
        for cutoff, bed_file in inferred_tract_files:
            precision, recall = cal_accuracy_samples(input.true_tracts, bed_file)
            acc_file = bed_file.rsplit('.', 1)[0] + ".accuracy"
            with open(acc_file, 'w') as f:
                f.write(f"{wildcards.scenario}\t{wildcards.replicate}\t{cutoff}\t{precision}\t{recall}\n")

    
        all_segments.to_csv(output.prob_file, index=False)

rule hmmix_artemis_accuracy_aggregated:
    input:
        prob_file = lambda wc: [
            os.path.join(
                results_folder,
                "hmmix_artemis",
                wc.scenario,
                str(rep),
                f"combined_segments.txt"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ],
        true_tracts=lambda wc: [
            os.path.join(
                data_folder,
                wc.scenario,
                str(rep),
                f"{scenario_configs[wc.scenario]['output_prefix']}.{rep}.true.tracts.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ],
    params:
        inferred_tracts=lambda wc: [
            os.path.join(
                results_folder,
                "hmmix_artemis",
                wc.scenario,
                str(rep),
                f"{scenario_configs[wc.scenario]['output_prefix']}.hybrid.{wc.cutoff}.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ],
        seq_len=lambda wc: int(scenario_configs[wc.scenario]["seq_len"])
    output:
        out=os.path.join(
            results_folder,
            "hmmix_artemis",
            "{scenario}",
            "hmmix_artemis.aggregated.{cutoff}.txt"
        ),
        metrics=os.path.join(
            results_folder,
            "hmmix_artemis",
            "{scenario}",
            "hmmix_artemis.aggregated.{cutoff}.metrics"
        )
    run:
        precision, recall, TP, FP, FN, TN = cal_metrics_multiple(
            input.true_tracts,
            params.inferred_tracts,
            params.seq_len
        )

        f1 = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0.0

        os.makedirs(os.path.dirname(output.out), exist_ok=True)

        with open(output.out, "w") as f:
            f.write(f"{wildcards.scenario}\t{wildcards.cutoff}\t{precision}\t{recall}\n")

        with open(output.metrics, "w") as f:
            f.write(
                f"{wildcards.scenario}\t{wildcards.cutoff}\t"
                f"{TP}\t{FP}\t{FN}\t{TN}\t"
                f"{precision}\t{recall}\t{f1}\n"
            )


rule hmmix_artemis_accuracy_summary:
    input:
        expand(
            os.path.join(
                results_folder,
                "hmmix_artemis",
                "{scenario}",
                "hmmix_artemis.aggregated.{cutoff}.txt"
            ),
            scenario=scenario_list,
            cutoff=cutoff_list
        )
    output:
        os.path.join(
            results_folder,
            "hmmix_artemis",
            "hmmix_artemis_accuracy_aggregated.txt"
        )
    shell:
        """
        cat {input} | sed '1i scenario\tcutoff\tprecision\trecall' > {output}
        """


rule hmmix_artemis_metrics_summary:
    input:
        expand(
            os.path.join(
                results_folder,
                "hmmix_artemis",
                "{scenario}",
                "hmmix_artemis.aggregated.{cutoff}.metrics"
            ),
            scenario=scenario_list,
            cutoff=cutoff_list
        )
    output:
        os.path.join(
            results_folder,
            "hmmix_artemis",
            "hmmix_artemis_metrics_aggregated.txt"
        )
    shell:
        """
        cat {input} | sed '1i scenario\tcutoff\tTP\tFP\tFN\tTN\tprecision\trecall\tF1' > {output}
        """



rule hmmix_all_done:
    input:
        get_all_hybrid_outputs
    output:
        os.path.join(skov_output_dir, "{scenario}", ".hmmix_done")
    shell:
        "touch {output}"