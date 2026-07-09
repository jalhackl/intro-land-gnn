import os

configfile: "run_relate_all.yaml"

#-----------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("{pref}_{chrom}.{ext}", chrom=config["chrom"], pref=config["pref"], ext=["anc", "mut"]),
#-----------------------------------------------------------------------------------------------------------------------
rule run_relate_rule_all:
    input:
        haps = "../../input/final/{pref}_{chrom}.haps.gz",
        sample = "../../input/final/{pref}_{chrom}.sample.gz",
        annot = "../../input/final/{pref}_{chrom}.annot",
        r_map = lambda wc: f"{config['rmaps_paths'][wc.chrom]}"
    output:
        anc = "{pref}_{chrom}.anc",
        mut = "{pref}_{chrom}.mut",
    log:
        '../../logs/run_relatemodeall_{pref}_{chrom}.log'
    params:
        prefix = "{pref}_{chrom}",
        path_to_relate = config["path_to_relate"],
        mu = config["mu"],
        haplotype_n = config["haplotype_n"]
    shell:
        """
         {params.path_to_relate}/bin/Relate \
          --mode All \
          -m {params.mu} \
          -N {params.haplotype_n} \
          --haps {input.haps} \
          --sample {input.sample} \
          --map {input.r_map} \
          --annot {input.annot} \
          -o {params.prefix} 
        """
#-----------------------------------------------------------------------------------------------------------------------
