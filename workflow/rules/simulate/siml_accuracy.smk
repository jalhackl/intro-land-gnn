import os
sys.path.append("workflow/scripts")
from evaluate_utils import cal_accuracy, cal_accuracy_multiple, cal_metrics_multiple
#-----------------------------------------------------------------------------------------------------------------------

rule siml_hmmix_accuracy:
    input:
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
        prob_file = rules.run_hmmix.output.prob_file,
    output:
        tracts = output_dir + output_prefix + "/{cutoff}.bed",
    resources:
        partition = "basic",
        mem_gb = 32,
        cpus = 16,
    conda:
        "../../../envs/hmmix-env.yaml",
    log:
        "logs/hmmix/simlaccuracy.{seed}.log",
    benchmark:
        "benchmarks/hmmix/simlaccuracy.{seed}.benchmark.txt",
    params:
        skov_output_dir = skov_output_dir,
        ref_set = ref_set,
        ref_id = config["ref_id"],
        src_id = config["src_id"],
        output_prefix = config["output_prefix"],
        cutoff_num = config["cutoff_num"],
        demog_id = demog_id,
        nref = nref,
        ntgt = ntgt
    shell:
        """    
        python  workflow/rules/hmmix/hmmix_siml_accuracy.py \
            --in_true_tracts '{input.true_tracts}' \
            --skov_output_dir '{params.skov_output_dir}' \
            --seed {wildcards.seed} \
            --out_prob_file '{input.prob_file}' \
            --src_id '{params.src_id}' \
            --output_prefix '{params.output_prefix}' \
            --cutoff_num {params.cutoff_num} \
            --demog_id '{params.demog_id}' \
            --nref {params.nref} \
            --ntgt {params.ntgt}  
        """


#-----------------------------------------------------------------------------------------------------------------------

rule sstar_process_output:
    input:
        quantiles = rules.sstar_threshold.output.quantiles,
        inferred_tracts = rules.sstar_output_to_bed.inferred_tracts,
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
    output:
        accuracy = "results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}" + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out.accuracy",
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    log:
        "logs/sstar/sstarprocout.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.log",
    benchmark:
        "benchmarks/sstar/sstarprocout.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.benchmark.txt",
    run:
        precision, recall = cal_accuracy(input.true_tracts, input.inferred_tracts)
        with open(output.accuracy, 'w') as o:
            o.write(f'{wildcards.demog}\t{wildcards.scenario}\tnref_{wildcards.nref}_ntgt_{wildcards.ntgt}\t{wildcards.quantile}\t{precision}\t{recall}\n')

 
rule accuracy_summary:
    input:
        accuracy_files = expand(sstar_output_dir + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out.accuracy",seed=seed_list,scenario=scenario_list, quantile=quantile_list),
    output:
        accuracy_table = os.path.join(sstar_output_dir + "/sstar_1src_accuracy.txt"),
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    log:
        "logs/sstar/accuracy.log",
    benchmark:
        "benchmarks/sstar/accuracy.benchmark.txt",
    shell:
        """
        cat {input.accuracy_files} | sed '1idemography\\tscenario\\tsample\\tcutoff\\tprecision\\trecall' > {output.accuracy_table}
        """
#-----------------------------------------------------------------------------------------------------------------------

rule sprime_process_output:
    input:
        scores = rules.sprime_run.output.score,
        inferred_tracts = rules.sprime_output_to_bed.inferred_tracts,
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
    output:
        accuracy = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.accuracy"),
    log:
        "logs/sprime/sprimeprocout.{demog}.{nref}.{ntgt}.{seed}/{threshold}.log",
    benchmark:
        "benchmarks/sprime/sprimeprocout.{demog}.{nref}.{ntgt}.{seed}/{threshold}.benchmark.txt",
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    run:
       precision, recall = cal_accuracy(input.true_tracts, input.inferred_tracts)
       with open(output.accuracy, 'w') as o:
           o.write(f'{wildcards.demog}\tnref_{wildcards.nref}_ntgt_{wildcards.ntgt}\t{wildcards.threshold}\t{precision}\t{recall}\n')

rule sprime_accuracy:
    input:
        accuracy_files = expand(sprime_output_dir + "/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.accuracy", demog=demog_id, nref=nref, ntgt=ntgt, seed=seed_list, threshold=threshold_list),
    output:
        accuracy_table = os.path.join(sprime_output_dir + "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/sprime_accuracy.txt"),
    log:
        "logs/sprime/sprime_accuracy.{demog}.{nref}.{ntgt}.{seed}.log",
    benchmark:
        "benchmarks/sprime/sprime_accuracy.{demog}.{nref}.{ntgt}.{seed}.benchmark.txt",
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    shell:
        """
        cat {input.accuracy_files} | sed '1idemography\\tscenario\\tsample\\tcutoff\\tprecision\\trecall' > {output.accuracy_table}
        """
#-----------------------------------------------------------------------------------------------------------------------
