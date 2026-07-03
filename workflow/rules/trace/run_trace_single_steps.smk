import os
import sys

sys.path.append("workflow/scripts")

from evaluate_utils import cal_accuracy, cal_accuracy_multiple, cal_metrics_multiple




# ------------------------------------------------------------
# STEP 1: build index file (no subprocess, no $(...), no hacks)
# ------------------------------------------------------------
rule build_trace_indices:
    input:
        ref=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.ref.ind.list"
        ),
        tgt=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tgt.ind.list"
        )

    output:
        idx=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tracehmm.idx"
        )

    params:
        ploidy=lambda wc: scenario_configs[wc.scenario].get("ploidy", 2)

    conda:
        "../../envs/trace-env.yaml"

    log:
        os.path.join("logs", "{scenario}", "{rep}", "{prefix}.{rep}.build_idx.log")

    shell:
        """
        python workflow/scripts/resolve_indices_relate.py \
            {input.ref} {input.tgt} {params.ploidy} \
            > {output.idx} 2> {log}
        """


# ------------------------------------------------------------
# STEP 2: run tracehmm extraction
# ------------------------------------------------------------
rule tracehmm_extract:
    input:
        trees=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.relate.trees"
        ),
        idx=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tracehmm.idx"
        )

    output:
        extracted=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{rep}",
            "{prefix}.{rep}.tracehmm.extract"
        )

    params:
        t=20000

    conda:
        "../../envs/trace-env.yaml"

    log:
        os.path.join("logs", "{scenario}", "{rep}", "{prefix}.{rep}.tracehmm_extract.log")

    shell:
        """
        python -m tracehmm.extract_cli \
            -f {input.trees} \
            -t {params.t} \
            -i {input.idx} \
            -o {output.extracted} \
            > {log} 2>&1
        """