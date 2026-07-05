from run_trace_pipeline import main

main(
    trees=snakemake.input.trees,
    ref_list=snakemake.input.ref_list,
    tgt_list=snakemake.input.tgt_list,
    scenario=snakemake.wildcards.scenario,
    rep=snakemake.wildcards.rep,
    prefix=snakemake.wildcards.prefix,
    out_base=snakemake.params.out_base,
    curr_folder=snakemake.params.curr_folder,
    t=snakemake.params.t,
    chrom=snakemake.params.chrom,
    posterior_threshold=snakemake.params.thresholds,
    #posterior_threshold=snakemake.wildcards.threshold,
    physical_length_threshold=snakemake.params.physical_length_threshold,
    genetic_distance_threshold=snakemake.params.genetic_distance_threshold,
)


open(snakemake.output.done, "w").write("done\n")