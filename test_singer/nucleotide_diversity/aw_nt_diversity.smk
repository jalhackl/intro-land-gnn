import os

typ=['wolf_n56','n82.nooutgroup']
chroms=['chr1',  'chr2',  'chr3',  'chr4',  'chr5',  'chr6',  'chr7',  'chr8',  'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18','chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27','chr28', 'chr29', 'chr30', 'chr31', 'chr32', 'chr33', 'chr34', 'chr35', 'chr36','chr37', 'chr38']
#-----------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz.tbi",typ=typ),
        expand("genetic_diversity/nucl_diversity/{typ}_{chrom}_pi",typ=typ, chrom=chroms),
#-----------------------------------------------------------------------------------------------------------------------
rule tabix_polarised:
    input:
        'anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz',
    output:
        'anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz.tbi',
    threads: 10
    shell:
        """
        tabix -s1 -b2 -e2 {input}
        """
#-----------------------------------------------------------------------------------------------------------------------
rule calc_nucl_diversity:
    input:
        'anc_allele/flip_alleles/aw_infosc0.8_maf0.01.autosomes.{typ}_aa_polarised.vcf.gz',
    output:
        'genetic_diversity/nucl_diversity/{typ}_{chrom}_pi'
    params:
        script_dir = "scripts/downstream/final_dataset/arg",
    shell:
        """
        python {params.script_dir}/aw_nt_diversity.py \
        --input_vcf {input} \
        --chrom {wildcards.chrom} \
        --pref {wildcards.typ}
        """
#-----------------------------------------------------------------------------------------------------------------------
