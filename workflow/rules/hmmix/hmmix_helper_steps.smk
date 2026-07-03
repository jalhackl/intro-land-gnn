import os
#---------------------------------------------------------------------------------------------------------------------
# Plot coalescent distribution to the outgroup for introgressed and ingroup states

rule hmmix_coalesc_dist:
    output:
        skov_helper_output_dir + "/{yaml}.pdf",
    resources:
        partition = "basic",
        mem_gb = 32,
        cpus = 6,
    conda:
        "../../../envs/hmmix-env.yaml",
    log:
        "logs/hmmix_helper/simulate_demog.{yaml}.log",
    params:
        yaml = config["yaml"],
        outgroupsize = config.get("outgroupsize", 500),
    shell:
        """    
        python resources/skov_helper_scripts/SimulateDemography/sim.py \
        -demography={params.yaml} \
        -iterations=10000 \
        -outgroupsize={params.outgroupsize}
        """
#---------------------------------------------------------------------------------------------------------------------
# default 10kb chr, sampled 1000 times -> 10Mb genome
rule hmmix_outgroup_size_explore:
    output:
        skov_helper_output_dir + "/outgroup_size_{yaml}.pdf",
    resources:
        partition = "basic",
        mem_gb = 32,
        cpus = 16,
    conda:
        "../../../envs/hmmix-env.yaml",
    log:
        "logs/hmmix_helper/outgroup_size_explore.{seed}.log",
    params:
        yaml = config["yaml"],
        outgroupsize = config.get("outgroupsize", 50),
        chrom_size = config.get("chrom_size", 10000),
        iterations = config.get("iterations", 1000),
    shell:
        """    
        python resources/skov_helper_scripts/Determine_outgroup_size/outgroup_efficientcy_github.py \
         -outgroup_size={params.outgroupsize} \
         -chrom_size={params.chrom_size} \
         -iterations={params.iterations} \
         -demography={params.yaml}.yaml \
         -outfile={output}
        """
#---------------------------------------------------------------------------------------------------------------------
