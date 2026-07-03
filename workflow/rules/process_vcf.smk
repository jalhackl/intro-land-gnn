checkpoint get_biallelic_vcf:
    input:
        vcf=lambda wc: os.path.join(
            main_config["test_data_folder"],
            wc.scenario,
            str(wc.replicate),
            f"{wc.output_prefix}.vcf"
        )
    output:
        vcf=os.path.join(
            main_config["test_data_folder"],
            "{scenario}",
            "{replicate}",
            "{output_prefix}.vcf.gz"
        )
    resources:
        time_min=120,
        mem_mb=5000,
        cpus=1
    shell:
        r"""
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {input.vcf}
        """
