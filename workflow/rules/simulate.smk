import os

rule simulate_gaishisim:
    output:
        # Stable completion marker
        done = touch(os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            ".simulate_done"
        ))

    resources:
        cpus=32,
        mem_gb=64

    log:
        "logs/{scenario}_simulate.log"

    params:
        threads=lambda wc, resources: resources.cpus * 2,

        output_dir=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario
        ),

        demes_file=lambda wc: scenario_configs[wc.scenario]['demes_file'],
        nref=lambda wc: scenario_configs[wc.scenario]['nref'],
        ntgt=lambda wc: scenario_configs[wc.scenario]['ntgt'],
        ref_id=lambda wc: scenario_configs[wc.scenario]['ref_id'],
        tgt_id=lambda wc: scenario_configs[wc.scenario]['tgt_id'],
        src_id=lambda wc: scenario_configs[wc.scenario]['src_id'],
        seq_len=lambda wc: scenario_configs[wc.scenario]['seq_len'],
        nrep=lambda wc: scenario_configs[wc.scenario]['nrep'],
        output_prefix=lambda wc: scenario_configs[wc.scenario]['output_prefix'],

        phased_flag=lambda wc:
            "--phased" if scenario_configs[wc.scenario].get("phased", False) else "",

        mutation_model=lambda wc:
            f"--mutation_model {scenario_configs[wc.scenario]['mutation_model']}"
            if "mutation_model" in scenario_configs[wc.scenario] else "",

        mut_map=lambda wc:
            f"--mut_map {scenario_configs[wc.scenario]['mut_map']}"
            if "mut_map" in scenario_configs[wc.scenario] else "",

        rec_map=lambda wc:
            f"--rec_map {scenario_configs[wc.scenario]['rec_map']}"
            if "rec_map" in scenario_configs[wc.scenario] else "",

        L_mut_map=lambda wc:
            f"--L_mut_map {scenario_configs[wc.scenario]['L_mut_map']}"
            if "L_mut_map" in scenario_configs[wc.scenario] else "",

        cv_mut_map=lambda wc:
            f"--cv_mut_map {scenario_configs[wc.scenario]['cv_mut_map']}"
            if "cv_mut_map" in scenario_configs[wc.scenario] else "",

        batch_processing_flag=lambda wc:
            "--batch_processing" if main_config.get("batch_processing", False) else "",

        shuffle_rec_map_flag=lambda wc:
            "--shuffle_rec_map"
            if scenario_configs[wc.scenario].get("shuffle_rec_map", False) else "",

    shell:
        """
        gaishisim simulate \
            --demes {params.demes_file} \
            --nref {params.nref} \
            --ntgt {params.ntgt} \
            --ref-id {params.ref_id} \
            --tgt-id {params.tgt_id} \
            --src-id {params.src_id} \
            --seq-len {params.seq_len} \
            --replicate {params.nrep} \
            --nprocess {params.threads} \
            --output-prefix {params.output_prefix} \
            --output-dir {params.output_dir} \
            {params.phased_flag} \
            {params.mutation_model} \
            {params.mut_map} \
            {params.rec_map} \
            {params.L_mut_map} \
            {params.cv_mut_map} \
            {params.batch_processing_flag} \
            {params.shuffle_rec_map_flag} \
            > {log} 2>&1
        """