import os

chroms=['1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14', '15', '16', '17', '18','19', '20', '21', '22', '23', '24', '25', '26', '27','28', '29', '30', '31', '32', '33', '34', '35', '36','37', '38']
#---------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("args/singer/input/{pref}_{chrom}.{ext}", chrom=chroms, pref='aw_n82', ext=['vcf.gz', 'vcf']),
#-----------------------------------------------------------------------------------------------------------------------
# split per chrom
rule split_vcf_per_chrom:
    input:
        'anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.n82.nooutgroup_aa_polarised.vcf.gz',
    output:
        'args/singer/input/{pref}_{chrom}.vcf.gz',
    log:
        'args/singer/logs/split_vcf_per_chrom_{pref}_{chrom}.log'
    params:
        region = "chr{chrom}",
    threads: 10
    shell:
        """
        bcftools view -r {params.region} {input} -Oz -o {output} --threads {threads}        
        """
#---------------------------------------------------------------------------------------------------------------------
# uncompress vcf
rule uncompress_vcf:
    input:
        rules.split_vcf_per_chrom.output
    output:
        'args/singer/input/{pref}_{chrom}.vcf',
    log:
        'args/singer/logs/uncompress_{pref}_{chrom}.log'
    threads: 10
    shell:
        """
        bgzip -dc {input} > {output}       
        """
#---------------------------------------------------------------------------------------------------------------------
