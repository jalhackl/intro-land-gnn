import os

configfile: "mk_relate_input.yaml"

#-----------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("args/relate/input/orig/{pref}_{chrom}.vcf.gz", chrom=config["chrom"], pref=config["pref"]),
        expand("args/relate/input/haps/{pref}_{chrom}.{ext}", chrom=config["chrom"], pref=config["pref"], ext=["haps", "sample"]),
        expand("args/relate/input/final/{pref}_{chrom}.{ext1}", chrom=config["chrom"], pref=config["pref"], ext1=['haps.gz', 'sample.gz', "annot"]),
#-----------------------------------------------------------------------------------------------------------------------
# generated outside of smk
#rule outgroup_to_fasta:
#    input:
#        in_bcf = "impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.n82.bcf",
#        fasta = 'gts/filter_ref_fasta/canFam31_autosomes.fasta',
#    output:
#        'args/relate/input/cafam31_anc.fa'
#    params:
#        outgroup = config["outgroup"]
#    log:
#        'args/relate/logs/outgroup_to_fasta.log'
#    threads: 10
#    shell:
#        """
#        bcftools consensus  -f {input.fasta} -s {params.outgroup} {input.in_bcf} > {output}
#        """
#-----------------------------------------------------------------------------------------------------------------------
# split per chrom
rule split_vcf_per_chrom:
    input:
        'impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.n82.nooutgroup.vcf.gz',
    output:
        'args/relate/input/orig/{pref}_{chrom}.vcf.gz',
    log:
        'args/relate/logs/split_vcf_per_chrom_{pref}_{chrom}.log'
    params:
        region = "chr{chrom}",
    threads: 10
    shell:
        """
        bcftools view -r {params.region} {input} -Oz -o {output} --threads {threads}
        tabix -f {output}       
        """
#---------------------------------------------------------------------------------------------------------------------
# convert vcf -> haps/sample

rule vcf_to_haps:
    input:
        'args/relate/input/orig/{pref}_{chrom}.vcf.gz',
    output:
        haps = 'args/relate/input/haps/{pref}_{chrom}.haps',
        sample = 'args/relate/input/haps/{pref}_{chrom}.sample',
    log:
        'args/relate/logs/vcf_to_haps_{pref}_{chrom}.log'
    params:
        prefix = "args/relate/input/orig/{pref}_{chrom}",
        path_to_relate = config["path_to_relate"]
    shell:
        """
         {params.path_to_relate}/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps {output.haps} \
                 --sample {output.sample} \
                 -i {params.prefix} \
                --chr {wildcards.chrom}       
        """
#-----------------------------------------------------------------------------------------------------------------------
# without mask
rule mk_relate_input:
    input:
        haps = rules.vcf_to_haps.output.haps,
        sample = rules.vcf_to_haps.output.sample,
        anc = "args/relate/input/cafam31_anc.fa",
        poplabels = config["poplabels"]
    output:
        hapsgz = "args/relate/input/final/{pref}_{chrom}.haps.gz",
        samplegz = "args/relate/input/final/{pref}_{chrom}.sample.gz",
        annot = "args/relate/input/final/{pref}_{chrom}.annot",
    log:
        'args/relate/logs/mkrelateinput_{pref}_{chrom}.log'
    params:
        prefix = "args/relate/input/final/{pref}_{chrom}",
        path_to_relate = config["path_to_relate"]
    shell:
        """
         {params.path_to_relate}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
                 --haps {input.haps} \
                 --sample {input.sample} \
                 --ancestor {input.anc} \
                 --poplabels {input.poplabels} \
                 -o {params.prefix}      
        """
#-----------------------------------------------------------------------------------------------------------------------