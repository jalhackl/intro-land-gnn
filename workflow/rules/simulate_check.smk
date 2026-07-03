import os

checkpoint simulate_gaishisim:
    output:
        done=os.path.join(main_config["test_data_folder"], "{scenario}", ".simulate_done")
    resources:
        cpus=32,
        mem_gb=64
    log:
        "logs/{scenario}_simulate.log"
    params:
        threads_used=lambda wildcards, resources: resources.cpus * 2,
        output_dir=lambda wildcards: os.path.join(
            main_config["test_data_folder"],
            wildcards.scenario
        ),
        demes_file=lambda wildcards: scenario_configs[wildcards.scenario]['demes_file'],
        nref=lambda wildcards: scenario_configs[wildcards.scenario]['nref'],
        ntgt=lambda wildcards: scenario_configs[wildcards.scenario]['ntgt'],
        ref_id=lambda wildcards: scenario_configs[wildcards.scenario]['ref_id'],
        tgt_id=lambda wildcards: scenario_configs[wildcards.scenario]['tgt_id'],
        src_id=lambda wildcards: scenario_configs[wildcards.scenario]['src_id'],
        seq_len=lambda wildcards: scenario_configs[wildcards.scenario]['seq_len'],
        nrep=lambda wildcards: scenario_configs[wildcards.scenario]['nrep'],
        output_prefix=lambda wildcards: scenario_configs[wildcards.scenario]['output_prefix'],
        phased_flag=lambda wildcards: "--phased" if scenario_configs[wildcards.scenario].get("phased", False) else "",
        mutation_model=lambda wildcards: f"--mutation_model {scenario_configs[wildcards.scenario]['mutation_model']}" 
                                    if "mutation_model" in scenario_configs[wildcards.scenario] else "",
        mut_map=lambda wildcards: f"--mut_map {scenario_configs[wildcards.scenario]['mut_map']}" 
                                    if "mut_map" in scenario_configs[wildcards.scenario] else "",
        rec_map=lambda wildcards: f"--rec_map {scenario_configs[wildcards.scenario]['rec_map']}" 
                                    if "rec_map" in scenario_configs[wildcards.scenario] else "",
        L_mut_map=lambda wildcards: f"--L_mut_map {scenario_configs[wildcards.scenario]['L_mut_map']}"
                                    if "L_mut_map" in scenario_configs[wildcards.scenario] else "",
        cv_mut_map=lambda wildcards: f"--cv_mut_map {scenario_configs[wildcards.scenario]['cv_mut_map']}"
                                    if "cv_mut_map" in scenario_configs[wildcards.scenario] else "",
        batch_processing_flag=lambda wildcards: "--batch_processing" 
                                    if main_config.get("batch_processing", False) else "",
        shuffle_rec_map_flag=lambda wildcards: "--shuffle_rec_map" 
                                    if scenario_configs[wildcards.scenario].get("shuffle_rec_map", False) else "",
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
            --nprocess {params.threads_used} \
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

        touch {output.done}
        """
