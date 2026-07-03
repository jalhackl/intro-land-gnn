import os


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
            "{prefix}.{rep}.tracehmm.extract.summary.full.summary.txt"
        )

    params:
        out_base=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            wc.rep,
            f"{wc.prefix}.{wc.rep}.tracehmm.extract"
        ),
        t=20000,
        curr_folder = main_config["test_data_folder"]

    conda:
        "../../envs/trace-env.yaml"

    log:
        os.path.join(
            "logs",
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tracehmm_pipeline.log"
        )

    shell:
        """
        python workflow/rules/trace/run_trace_pipeline.py \
            {input.trees} \
            {input.ref_list} \
            {input.tgt_list} \
            {wildcards.scenario} \
            {wildcards.rep} \
            {wildcards.prefix} \
            {params.out_base} \
            {params.curr_folder} \
            {params.t} \
            > {log} 2>&1
        """