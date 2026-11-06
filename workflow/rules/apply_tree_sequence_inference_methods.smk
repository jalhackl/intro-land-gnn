import os

from gaishisim.utils.simulators.InferTskit import InferTskit
from gaishisim.utils.simulators.InferRelate import InferRelate
from gaishisim.utils.simulators.InferSinger import InferSinger


# ============================================================
# tsinfer
# ============================================================

rule apply_tsinfer:
    input:
        vcf=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.vcf"
        ),
        ts=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.ts"
        )

    output:
        tsinfer=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}_tsinfer.ts"
        )

    params:
        output_dir=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep)
        ),
        mutation_rate=lambda wc: scenario_configs[wc.scenario].get(
            "mut_rate",
            1e-8
        )

    resources:
        cpus=32,
        mem_gb=64

    log:
        os.path.join(
            "logs",
            "{scenario}",
            "{rep}",
            "{prefix}",
            "tsinfer.log"
        )

    run:
        os.makedirs(params.output_dir, exist_ok=True)

        infer = InferTskit(
            vcf_file=input.vcf,
            output_prefix=f"{wildcards.prefix}.{wildcards.rep}",
            output_dir=params.output_dir,

            infer_trees=True,
            remove_zarr=True,

            ancestral_mode="from_ts",
            ts_file=input.ts,
            check_ancestral_state=True,

            apply_tsdate=True,
            tsdate_mutation_rate=params.mutation_rate,

            worker_processes=resources.cpus,
        )

        infer.run()

        if not os.path.exists(output.tsinfer):
            raise RuntimeError(
                f"tsinfer did not create expected output: "
                f"{output.tsinfer}"
            )


# ============================================================
# Relate
# ============================================================

rule apply_relate:
    input:
        vcf=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.vcf"
        ),
        ts=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.ts"
        )

    output:
        trees=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.relate.trees"
        ),
        # additional output for metadata and reestimation
                trees_meta=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.relate_reestimated.metadata.trees"
        ),
                        trees_reest_meta=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.relate_reestimated.trees"
        )

    params:
        output_prefix=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.relate"
        ),

        demes_file=lambda wc: scenario_configs[wc.scenario]["demes_file"],

        seq_len=lambda wc: scenario_configs[wc.scenario]["seq_len"],

        mutation_rate=lambda wc: scenario_configs[wc.scenario].get(
            "mut_rate",
            1e-8
        ),

        recombination_rate=lambda wc: scenario_configs[wc.scenario].get(
            "rec_rate",
            1e-8
        ),

        relate_memory=lambda wc: scenario_configs[wc.scenario].get(
            "relate_memory",
            5
        ),

        relate_mode=lambda wc: scenario_configs[wc.scenario].get(
            "relate_mode",
            "parallel"
        )

    resources:
        cpus=32,
        mem_gb=64

    log:
        os.path.join(
            "logs",
            "{scenario}",
            "{rep}",
            "{prefix}",
            "relate.log"
        )

    run:
        os.makedirs(
            os.path.dirname(params.output_prefix),
            exist_ok=True
        )

        infer = InferRelate(
            vcf_file=input.vcf,
            output_prefix=params.output_prefix,

            recombination_rate=params.recombination_rate,
            sequence_length=params.seq_len,
            mutation_rate=params.mutation_rate,

            demes_file=params.demes_file,

            Ne=None,
            Ne_estimate_method="diversity",

            mode=params.relate_mode,
            threads=resources.cpus,
            memory=params.relate_memory,

            original_ts_file=input.ts,
        )

        infer.run()

        if not os.path.exists(output.trees):
            raise RuntimeError(
                f"Relate did not create expected output: "
                f"{output.trees}"
            )


# ============================================================
# SINGER
# ============================================================


rule apply_singer:
    input:
        vcf=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.vcf"
        ),
        ts=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.ts"
        ),
        poplabels=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.poplabels"
        )

    output:
        trees=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.singer.trees"
        )

    params:
        output_prefix=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.rep),
            f"{wc.prefix}.{wc.rep}.singer"
        ),
        mutation_rate=lambda wc: scenario_configs[wc.scenario].get(
            "mut_rate", 1e-8
        ),
        seq_len=lambda wc: scenario_configs[wc.scenario]["seq_len"],
        rec_map=lambda wc: scenario_configs[wc.scenario].get(
            "rec_map", None
        ),
        mut_map=lambda wc: scenario_configs[wc.scenario].get(
            "mut_map", None
        )

    resources:
        cpus=8,
        mem_gb=64,
		time=1400

    log:
        os.path.join(
            "logs",
            "{scenario}",
            "{rep}",
            "{prefix}",
            "singer.log"
        )

    run:
        os.makedirs(
            os.path.dirname(params.output_prefix),
            exist_ok=True
        )

        infer = InferSinger(
            vcf_file=input.vcf,
            output_prefix=params.output_prefix,
            mutation_rate=params.mutation_rate,
            sequence_length=params.seq_len,
            recombination_map=params.rec_map,
            mutation_map=params.mut_map,
            original_ts=input.ts,
            convert_to_trees=True,
            Ne=None,
            Ne_estimate_method="diversity",

            # Populations/sample metadata from the simulator
            poplabels_file=input.poplabels,
            add_metadata=True,
        )

        infer.run()

        # InferSinger creates {output_prefix}_0.trees
        singer_tree = f"{params.output_prefix}_0.trees"

        if not os.path.exists(singer_tree):
            raise RuntimeError(
                f"SINGER did not create expected intermediate output: "
                f"{singer_tree}"
            )

        # Rename it according to Snakemake output
        os.replace(singer_tree, output.trees)

        if not os.path.exists(output.trees):
            raise RuntimeError(
                f"SINGER did not create expected output: {output.trees}"
            )
