#-----------------------------------------------------------------------------------------------------------------------
#rule all:
#    input:
#        expand("resources/{hgpref}/{hgpref}_{ext}", hgpref=['hg19', 'hg38'], ext=['ancestral.tar.gz', 'refgenome.tar.gz', 'strick_callability_mask.bed']),
#        directory("resources/skov_helper_scripts/Determine_outgroup_size"),
#        directory("resources/skov_helper_scripts/SimulateDemography")
#-----------------------------------------------------------------------------------------------------------------------
# download skov files - https://github.com/LauritsSkov/Introgression-detection?tab=readme-ov-file
# anc genomes, callability masks, ref genomes
# curl  https://zenodo.org/api/records/11212339  | grep 'links'
rule get_hg19_hg38_files:
    log:
        "logs/get_{hgpref}files.log"
    output:
        "resources/{hgpref}/{hgpref}_ancestral.tar.gz",
        "resources/{hgpref}/{hgpref}_refgenome.tar.gz",
        "resources/{hgpref}/{hgpref}_strick_callability_mask.bed",
    params:
        outdir = "resources/{hgpref}"
    resources: nodes=1, ntasks=1, time_min=240, mem_gb=100, cpus=4
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        wget https://zenodo.org/api/records/11212339/files/{wildcards.hgpref}_ancestral.tar.gz
        wget https://zenodo.org/api/records/11212339/files/{wildcards.hgpref}_refgenome.tar.gz
        wget https://zenodo.org/api/records/11212339/files/{wildcards.hgpref}_strick_callability_mask.bed
        """
#-----------------------------------------------------------------------------------------------------------------------
rule get_hmmix_helper_scripts:
    output:
        directory("resources/skov_helper_scripts/Determine_outgroup_size"),
    resources: nodes=1, ntasks=1, time_min=60, mem_gb=30, cpus=1
    log:
        "logs/gitclone_skovscripts.log"
    params:
        outdir = "resources/skov_helper_scripts"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        git clone https://github.com/LauritsSkov/Determine_outgroup_size
        """
#-----------------------------------------------------------------------------------------------------------------------
rule get_hmmix_helper_scripts1:
    output:
        directory("resources/skov_helper_scripts/SimulateDemography")
    resources: nodes=1, ntasks=1, time_min=60, mem_gb=30, cpus=1
    log:
        "logs/gitclone_skovscripts1.log"
    params:
        outdir = "resources/skov_helper_scripts"
    shell:
        """
        cd {params.outdir}
        git clone https://github.com/LauritsSkov/SimulateDemography
        """
#-----------------------------------------------------------------------------------------------------------------------



