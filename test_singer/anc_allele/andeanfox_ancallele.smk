import os

typ=['wolf_n56','n82.nooutgroup']
#-----------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("anc_allele/andeanfox_anc.gz"),
        expand("anc_allele/aa.hdr.txt"),
        expand("anc_allele/aw_infosc0.8_maf0.01.autosomes.{typ}_aa.vcf.gz",typ=typ),
        expand("anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz",typ=typ),
#-----------------------------------------------------------------------------------------------------------------------
rule get_outgroup_alleles:
    output:
        'anc_allele/andeanfox_anc',
    params:
        Rscript_dir = "scripts/downstream/final_dataset/arg",
    shell:
        """
        Rscript --vanilla {params.Rscript_dir}/andeanfox_ancallele.R
        """
#-----------------------------------------------------------------------------------------------------------------------
rule tabix_outgroup_alleles:
    input:
        in_bcf = rules.get_outgroup_alleles.output,
    output:
        'anc_allele/andeanfox_anc.gz'
    threads: 10
    shell:
        """
        bgzip {input}
        tabix -s1 -b2 -e2 {output}
        """
#-----------------------------------------------------------------------------------------------------------------------
rule make_aa_header:
    output:
        'anc_allele/aa.hdr.txt'
    threads: 10
    shell:
        """
        echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > {output}
        """
#-----------------------------------------------------------------------------------------------------------------------
rule annotate_vcf_aa:
    input:
        in_bcf = "impute/joint_impute/2905g_modern/filt/CanFam31/concat/norel/aw_infosc0.8_maf0.01.autosomes.{typ}.bcf",
        aa = 'anc_allele/andeanfox_anc.gz',
        hdr = 'anc_allele/aa.hdr.txt'
    output:
        'anc_allele/aw_infosc0.8_maf0.01.autosomes.{typ}_aa.vcf.gz'
    threads: 10
    shell:
        """
        bcftools annotate \
        -a {input.aa} \
        -h {input.hdr} \
        -c CHROM,POS,REF,ALT,INFO/AA \
        {input.in_bcf} \
        -Oz -o {output}
        bcftools index -f {output}
        """
#-----------------------------------------------------------------------------------------------------------------------
rule flip_alleles_vcf_aa:
    input:
        'anc_allele/aw_infosc0.8_maf0.01.autosomes.{typ}_aa.vcf.gz',
    output:
        'anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz'
    params:
        script_dir = "scripts/downstream/final_dataset/arg",
    shell:
        """
        python {params.script_dir}/flip_alleles.py \
        --input_vcf {input} \
        --output_vcf {output}
        bcftools index -f {output}
        """
#-----------------------------------------------------------------------------------------------------------------------