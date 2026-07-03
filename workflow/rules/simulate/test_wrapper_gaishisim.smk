import os
configfile: "config/simulate/test_gaishisim_config.yaml"

KEYS = list(config["scensimulate"].keys())
print(KEYS)
nrep = config["nrep"]
#rep = range(nrep)
#print(rep)

rule all:
    input:
        expand("sim/{sim_output_prefix}/{rep}/{sim_output_prefix}.{ext}", sim_output_prefix=KEYS, rep=range(nrep), ext=['ts', 'vcf.gz', 'ref.ind.list', 'tgt.ind.list'])


include: "test_gaishisim1.smk"

