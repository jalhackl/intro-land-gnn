# takes dictionary of the scenarios with the required params

rule simulate_gaishsim_1:
    output:
        ts = lambda wildcards: expand("sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.ts", sim_output_prefix=[wildcards.sim_output_prefix], rep=range(config["nrep"])),
        vcf = lambda wildcards: expand("sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.vcf", sim_output_prefix=[wildcards.sim_output_prefix], rep=range(config["nrep"])),
        ref = lambda wildcards: expand("sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.ref.ind.list", sim_output_prefix=[wildcards.sim_output_prefix], rep=range(config["nrep"])),
        tgt = lambda wildcards: expand("sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.tgt.ind.list", sim_output_prefix=[wildcards.sim_output_prefix], rep=range(config["nrep"])),
    log:
        "logs/simulate_test_data/sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.log",
    resources:
        cpus=32,
        mem_gb=64
    params:
        sim_output_prefix = "{sim_output_prefix}",
        sim_outputdir = "sim/{sim_output_prefix}",
        ploidy = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("ploidy", 2),
        seq_len = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("seq_len", 25000000),
        nrep = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("nrep", 10),
        mut_rate = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("mut_rate", 1.25e-8),
        rec_rate = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("rec_rate", 1.0e-8),
        ref_id = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix]["ref_id"],
        tgt_id = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix]["tgt_id"],
        src_id = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix]["src_id"],
        nref = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("nref", 50),
        ntgt = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix].get("ntgt", 10),
        demes_file = lambda wildcards: config["scensimulate"][wildcards.sim_output_prefix]["demes_file"],
        threads_used = lambda wildcards, resources: resources.cpus * 2,
        phased_flag=lambda wildcards: "--phased" if config["scensimulate"][wildcards.sim_output_prefix].get("phased", False) else "",
        mutation_model=lambda wildcards: f"--mutation_model {config['scensimulate'][wildcards.sim_output_prefix]['mutation_model']}" 
                                    if "mutation_model" in config["scensimulate"][wildcards.sim_output_prefix] else "",
        mut_map=lambda wildcards: f"--mut_map {config['scensimulate'][wildcards.sim_output_prefix]['mut_map']}" 
                                    if "mut_map" in config["scensimulate"][wildcards.sim_output_prefix] else "",
        rec_map=lambda wildcards: f"--rec_map {config['scensimulate'][wildcards.sim_output_prefix]['rec_map']}" 
                                    if "rec_map" in config["scensimulate"][wildcards.sim_output_prefix] else "",
        L_mut_map=lambda wildcards: f"--L_mut_map {config['scensimulate'][wildcards.sim_output_prefix]['L_mut_map']}"
                                    if "L_mut_map" in config["scensimulate"][wildcards.sim_output_prefix] else "",
        cv_mut_map=lambda wildcards: f"--cv_mut_map {config['scensimulate'][wildcards.sim_output_prefix]['cv_mut_map']}"
                                    if "cv_mut_map" in config["scensimulate"][wildcards.sim_output_prefix] else "",
        batch_processing_flag=lambda wildcards: "--batch_processing" 
                                    if main_config.get("batch_processing", False) else "",
        shuffle_rec_map_flag=lambda wildcards: "--shuffle_rec_map" 
                                    if config["scensimulate"][wildcards.sim_output_prefix].get("shuffle_rec_map", False) else "",
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
            --output-prefix {params.sim_output_prefix} \
            --output-dir {params.sim_output_dir} \
            {params.phased_flag} \
            {params.mutation_model} \
            {params.mut_map} \
            {params.rec_map} \
            {params.L_mut_map} \
            {params.cv_mut_map} \
            {params.batch_processing_flag} \
            {params.shuffle_rec_map_flag} \
        """

rule compress_sim_vcf:
    input:
        vcf = "sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.vcf",
    output:
        vcf = "sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.vcf.gz",
    resources:
        time_min=120,
        mem_mb=5000,
        cpus=1
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {input.vcf}
        """
