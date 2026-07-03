import os
import demes
import numpy as np
import pandas as pd
import sys

from sstar_additional_functions import *

sys.path.append("worflow/rules")
sys.path.append("worflow/rules/sstar")
from evaluate_utils import *


'''
rule all:
    input:
        sstar_output_dir + "/sstar_1src_accuracy.txt",
'''

data_folder = main_config["test_data_folder"]
results_folder = main_config["results_folder"]
scenario_list = list(scenario_configs.keys())
quantile_list   = sstar_config["quantile_list"]

# Required for simulation quantile summaries
snp_num_list    = np.arange(25, 705, 5)

sstar_output_dir_simulation = os.path.join(results_folder, "simulation")
sstar_output_dir = os.path.join(results_folder, "sstar")


# get replicates for a scenario after checkpoint finishes
def get_replicates_for_scenario(scenario):
    ck = checkpoints.simulate_gaishisim.get(scenario=scenario)
    scenario_dir = os.path.join(data_folder, scenario)
    replicates = sorted([
        d for d in os.listdir(scenario_dir)
        if os.path.isdir(os.path.join(scenario_dir, d)) and not d.startswith(".")
    ])
    return replicates

rule mut_rec_combination:
    output:
        rates = sstar_output_dir_simulation + "/{scenario}/snps/rates.combination",
    params:
        seq_len = sstar_config.get("ms_seq_len", 50000),
        mut_rate = lambda wc: float(scenario_configs[wc.scenario]['mut_rate']),
        rec_rate = lambda wc: float(scenario_configs[wc.scenario]['rec_rate']),
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    conda:
        "../../envs/sstar-env.yaml",
    log:
        "logs/sstar/mutrec.{scenario}.log",
    benchmark:
        "benchmarks/sstar/mutrec.{scenario}.benchmark.txt",
    shell:
        """    
        python workflow/rules/sstar/mut_rec.py \
         --seqlen {params.seq_len}  --mutrate {params.mut_rate} \
         --recrate {params.rec_rate} --outfile '{output.rates}'
        """

rule simulate_glm_data:
    input:
        rates = rules.mut_rec_combination.output.rates,
    output:
        ms = sstar_output_dir_simulation + "/{scenario}/snps/{snp_num}/sim1src.ms",
    log:
        "logs/sstar/msglm.{scenario}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/msglm.{scenario}.{snp_num}.benchmark.txt",
    params:
        nreps = sstar_config.get("ms_nreps", 20000),
        seq_len = sstar_config.get("ms_seq_len", 50000),
        ms_exec = sstar_config["ms_exec"],
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        tgt_id = scenario_configs[wildcards.scenario]['tgt_id']
        ref_id = scenario_configs[wildcards.scenario]['ref_id']

        ntgt = int(scenario_configs[wildcards.scenario]['ntgt'])
        nref = int(scenario_configs[wildcards.scenario]['nref'])
        demes_file = scenario_configs[wildcards.scenario]["no_intro_demes_file"]

        new_graph = demes.load(demes_file)
        demes_length = len(new_graph._deme_map)
        ref_index = 0
        tgt_index = 0

        for i, deme in enumerate(new_graph._deme_map):
            if deme == ref_id:
                ref_index = i
            if deme == tgt_id:
                tgt_index = i
        ms_samples = demes_length * [0]
        ms_samples[tgt_index] = 2
        ms_samples[ref_index] = nref * 2

        nsamp = 2*(1 + nref)

        ms_params = demes.to_ms(new_graph, N0=1000, samples=ms_samples)

        out_dir = os.path.dirname(output.ms)
        shell(
            "bash -c '"
            "set -e; "  # exit on error, but no pipefail
            "mkdir -p {out_dir} && "
            "cat {input.rates} | {params.ms_exec} {nsamp} {params.nreps} "
            "-t tbs -r tbs {params.seq_len} -s {wildcards.snp_num} {ms_params} > {output.ms}'"
        )

#-----------------------------------------------------------------------------------------------------------------------
# create ind file lists for ms simulations
rule ms_ind_lists:
    output:
        ss_ref = sstar_output_dir_simulation + "/{scenario}/sstarsim.ref.ind.list",
        ss_tgt = sstar_output_dir_simulation + "/{scenario}/sstarsim.tgt.ind.list",
    log:
        "logs/sstar/ms_ind_lists.{scenario}.log",
    benchmark:
        "benchmarks/sstar/ms_ind_lists.{scenario}.benchmark.txt",
    params:
        ploidy = lambda wc: int(scenario_configs[wc.scenario]['ploidy'])
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        ntgt = int(scenario_configs[wildcards.scenario]['ntgt'])
        nref = int(scenario_configs[wildcards.scenario]['nref'])

        nsamp = 2*(1 + nref)
        create_ind_lists(nsamp, output.ss_ref, output.ss_tgt, params.ploidy, ind_prefix="tsk_")

rule ms2vcf:
    input:
        ms = rules.simulate_glm_data.output.ms,
    output:
        vcf = sstar_output_dir_simulation + "/{scenario}/snps/{snp_num}/sim1src.vcf",
    log:
        "logs/sstar/ms2vcf.{scenario}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/ms2vcf.{scenario}.{snp_num}.benchmark.txt",
    params:
        seq_len = sstar_config.get("ms_seq_len", 50000),
        ploidy = lambda wc: int(scenario_configs[wc.scenario]['ploidy'])
    resources: 
        time = lambda wildcards: 360 if (int(wildcards.snp_num) < 340) else 1000,
        mem_gb = lambda wildcards: 10 if (int(wildcards.snp_num) < 340) else 100,
        cpus = 1
    threads: 1,
    run:
        ntgt = int(scenario_configs[wildcards.scenario]['ntgt'])
        nref = int(scenario_configs[wildcards.scenario]['nref'])

        nsamp = 2*(1 + nref)
        ms2vcf(input.ms, output.vcf, nsamp, params.seq_len, params.ploidy, ind_prefix="tsk_")
#-----------------------------------------------------------------------------------------------------------------------


rule cal_score:
    input:
        vcf = rules.ms2vcf.output.vcf,
	ref_list = rules.ms_ind_lists.output.ss_ref,
	tgt_list = rules.ms_ind_lists.output.ss_tgt,
    output:
        score = sstar_output_dir_simulation + "/{scenario}/snps/{snp_num}/sim1src.sstar.scores",
    params:
        seq_len = sstar_config.get("ms_seq_len", 50000),
    conda:
        "../../envs/sstar-env.yaml",
    benchmark:
        "benchmarks/sstar/sstarcalcscore.{scenario}.{snp_num}.benchmark.txt",
    log:
        "logs/sstar/sstarcalcscore.{scenario}.{snp_num}.log",
    resources:
        time = lambda wildcards: 360 if (int(wildcards.snp_num) < 340) else 1000,
        mem_gb = lambda wildcards: 10 if (int(wildcards.snp_num) < 340) else 100,
        cpus = 1
    threads: 1,
    shell:
        """
        sstar score --vcf {input.vcf} --ref {input.ref_list} --tgt {input.tgt_list} --output {output.score} --thread {threads} --win-len {params.seq_len} --win-step {params.seq_len}
        """


rule cal_quantile:
    input:
        score = rules.cal_score.output.score,
    output:
        quantile = sstar_output_dir_simulation + "/{scenario}/snps/{snp_num}/sim1src.sstar.quantile",
    params:
        sim_quantiles = quantile_list,
    log:
        "logs/sstar/sstarcalcquant.{scenario}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/sstarcalcquant.{scenario}.{snp_num}.benchmark.txt",
    resources: time_min=3000, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        df = pd.read_csv(input.score, sep="\t").dropna()
        if df.empty:
            # Write empty quantile file to avoid breaking downstream rules
            with open(output.quantile, 'w') as o:
                for q in params.sim_quantiles:
                    o.write(f'NA\t{wildcards.snp_num}\t{q}\n')
        else:
            mean_df = df.groupby(['chrom', 'start', 'end'], as_index=False)['S*_score'].mean().dropna()
            if mean_df['S*_score'].empty:
                with open(output.quantile, 'w') as o:
                    for q in params.sim_quantiles:
                        o.write(f'NA\t{wildcards.snp_num}\t{q}\n')
            #regular case
            else:
                scores = np.quantile(mean_df['S*_score'], params.sim_quantiles)
                with open(output.quantile, 'w') as o:
                    for i in range(len(scores)):
                        o.write(f'{scores[i]}\t{wildcards.snp_num}\t{params.sim_quantiles[i]}\n')


rule quantile_summary:
    input:
        res = expand(sstar_output_dir_simulation + "/{scenario}/snps/{snp_num}/sim1src.sstar.quantile", snp_num=snp_num_list, scenario=scenario_list),
    output:
        output_res = sstar_output_dir_simulation + "/{scenario}/quantile.1src.summary.txt",
    resources: time_min=3000, mem_mb=5000, cpus=1,
    threads: 1,
    log:
        "logs/sstar/sstarquantsum.{scenario}.log",
    benchmark:
        "benchmarks/sstar/sstarquantsum.{scenario}.benchmark.txt",
    shell:
        """
        cat {input.res} | sort -nk 2,2 | sed '1iS*_score\\tSNP_number\\tquantile\n' > {output.output_res}
        """


rule sstar_score:
    input:
        vcf=lambda wc: checkpoints.get_biallelic_vcf.get(
            scenario=wc.scenario,
            replicate=wc.replicate,
            output_prefix=f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}"
        ).output.vcf,
        ref_ind=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            str(wc.replicate),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.ref.ind.list"
        ),
        tgt_ind=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            str(wc.replicate),
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.tgt.ind.list"
        )
    output:
        score = sstar_output_dir_simulation + "/{scenario}/{replicate}/sstar.1src.out.score",
    params:
        seq_len = sstar_config.get("score_seq_len", 50000),
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    conda:
        "../../envs/sstar-env.yaml",
    log:
        "logs/sstar/sstarscore.{scenario}{replicate}.log",
    benchmark:
        "benchmarks/sstar/sstarscore.{scenario}{replicate}.benchmark.txt",
    shell:
        """
        sstar score --vcf {input.vcf} --ref {input.ref_ind} --tgt {input.tgt_ind} --output {output.score} --thread {threads} --win-len {params.seq_len} --win-step 10000
        """


rule sstar_threshold:
    input:
        score = rules.sstar_score.output.score,
        summary = sstar_output_dir_simulation + "/{scenario}/quantile.1src.summary.txt",
    output:
        quantiles = sstar_output_dir_simulation + "/{scenario}/{replicate}/sstar.1src.quantile.{quantile}.out",

        #quantiles = output_dir + "inference/sstar/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out",
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    conda:
        "../../envs/sstar-env.yaml",
    log:
        "logs/sstar/sstarthresh.{scenario}.{replicate}.{quantile}.log",
    benchmark:
        "benchmarks/sstar/sstarthresh.{scenario}.{replicate}.{quantile}.benchmark.txt",
    shell:
        """
        sstar threshold --score {input.score} --sim-data {input.summary} --quantile {wildcards.quantile} --output {output.quantiles}
        """


#Rule: Process sstar output to generate inferred tracts and accuracy
rule sstar_process_output:
    input:
        quantiles = rules.sstar_threshold.output.quantiles,
        true_tracts=lambda wc: os.path.join(
            data_folder,
            wc.scenario,
            wc.replicate,
            f"{scenario_configs[wc.scenario]['output_prefix']}.{wc.replicate}.true.tracts.bed"
        )
    output:
        inferred_tracts=os.path.join(results_folder, "sstar", "{scenario}", "{replicate}", "sstar.1src.quantile.{quantile}.out.bed"),
        accuracy=os.path.join(results_folder, "sstar", "{scenario}", "{replicate}", "sstar.1src.quantile.{quantile}.out.accuracy")
    threads: 1
    resources:
        time_min=60,
        mem_mb=2000,
        cpus=1
    run:
        process_sstar_1src_output(input.quantiles, output.inferred_tracts)
        precision, recall = cal_accuracy_samples(input.true_tracts, output.inferred_tracts)
        with open(output.accuracy, 'w') as f:
            f.write(f"{wildcards.scenario}\t{wildcards.replicate}\t{wildcards.quantile}\t{precision}\t{recall}\n")



#Function to get all accuracy files 
def get_all_accuracy_files(wildcards):
    all_files = []
    for scenario in scenario_list:
        replicates = get_replicates_for_scenario(scenario)
        for rep in replicates:
            for quantile in quantile_list:
                all_files.append(os.path.join(
                 results_folder, "sstar", scenario, rep, f"sstar.1src.quantile.{quantile}.out.accuracy"
                ))
    return all_files


#Rule: Aggregate all accuracy results into one summary table
rule sstar_accuracy:
    input:
        accuracy_files=get_all_accuracy_files
    output:
        accuracy_table=os.path.join(results_folder, "sstar", "sstar_accuracy.txt")
    shell:
        """
        cat {input.accuracy_files} | sed '1i scenario\treplicate\tcutoff\tprecision\trecall' > {output.accuracy_table}
        """

#trigger all sprime runs & processing
rule sstar_all:
    input:
        get_all_accuracy_files



rule sstar_accuracy_aggregated:
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
        #inferred_tracts = rules.sstar_process_output.output.inferred_tracts
        inferred_tracts=lambda wc: [
            os.path.join(
                results_folder, "sstar", wc.scenario, rep,
                f"sstar.1src.quantile.{wc.quantile}.out.bed"
            )
            for rep in get_replicates_for_scenario(wc.scenario)
        ]
    params:
        seq_len=lambda wc: int(scenario_configs[wc.scenario]["seq_len"])
    output:
        out=os.path.join(results_folder, "sstar", "{scenario}", "sstar.aggregated.{quantile}.txt"),
        metrics=os.path.join(results_folder, "sstar", "{scenario}", "sstar.aggregated.{quantile}.metrics")
    run:
        print("Scenario:", wildcards.scenario)
        print("Threshold:", wildcards.quantile)
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
                f"{wildcards.quantile}\t"
                f"{precision}\t"
                f"{recall}\n"
            )

        # Write detailed metrics including F1-score
        with open(output.metrics, "w") as f:
            f.write(
               # "scenario\tthreshold\tTP\tFP\tFN\tTN\tprecision\trecall\tF1\n"
               f"{wildcards.scenario}\t"
               f"{wildcards.quantile}\t"
               f"{total_true_positive}\t"
               f"{total_false_positive}\t"
               f"{total_false_negative}\t"
               f"{total_true_negative}\t"
               f"{precision}\t"
               f"{recall}\t"
               f"{f1}\n"
           )



rule sstar_accuracy_summary:
    input:
        expand(os.path.join(results_folder, "sstar", "{scenario}", "sstar.aggregated.{quantile}.txt"),
               scenario=scenario_list, quantile=quantile_list)
    output:
        os.path.join(results_folder, "sstar", "sstar_accuracy_aggregated.txt")
    shell:
        """
        cat {input} | sed '1i scenario\tthreshold\tprecision\trecall' \
            > {output}
        """



rule sstar_metrics_summary:
    input:
        expand(os.path.join(results_folder, "sstar", "{scenario}", "sstar.aggregated.{quantile}.metrics"),
               scenario=scenario_list, quantile=quantile_list)
    output:
        os.path.join(results_folder, "sstar", "sstar_metrics_aggregated.txt")
    shell:
        """
        cat {input} | sed '1i scenario\tthreshold\tTP\tFP\tFN\tTN\tF1' \
            > {output}
        """


