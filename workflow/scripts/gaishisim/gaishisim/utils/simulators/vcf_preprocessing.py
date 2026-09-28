import pandas as pd
from cyvcf2 import VCF, Writer
import os
import tempfile
import shutil


def _prepare_output(input_vcf, output_vcf):
    """
    Prepare a safe output path.

    If input and output are identical, create a temporary file with
    the same compression type as the requested output.

    Returns:
        write_path: path where cyvcf2 should write
        final_path: path that should contain the final result
        temporary: whether replacement is needed afterwards
    """

    same_file = (
        os.path.abspath(input_vcf) == os.path.abspath(output_vcf)
    )

    if same_file:

        output_vcf = os.path.abspath(output_vcf)

        if output_vcf.endswith(".vcf.gz"):
            suffix = ".vcf.gz"
        elif output_vcf.endswith(".vcf"):
            suffix = ".vcf"
        else:
            suffix = os.path.splitext(output_vcf)[1]

        tmp = tempfile.NamedTemporaryFile(
            suffix=suffix,
            dir=os.path.dirname(output_vcf),
            delete=False
        )
        tmp.close()

        return tmp.name, output_vcf, True

    return output_vcf, output_vcf, False



def _finalize_output(write_path, final_path, temporary):
    """
    Replace original file if temporary output was used.
    """

    if temporary:
        shutil.move(write_path, final_path)




def flip_to_ancestral_ref(input_vcf, output_vcf):
    """
    Flip REF/ALT alleles so that AA is always REF.

    Requires an AA INFO field.

    Parameters
    ----------
    input_vcf : str
        Input VCF(.vcf.gz).
    output_vcf : str
        Output VCF(.vcf.gz).
    """

    write_path, final_path, temporary = _prepare_output(
        input_vcf, output_vcf
    )

    vcf = VCF(input_vcf)
    writer = Writer(write_path, vcf)

    for record in vcf:

        ancestral = record.INFO.get("AA")

        if (
            ancestral not in (None, ".", "")
            and len(record.ALT) == 1
            and record.ALT[0] == ancestral
        ):

            # Swap REF and ALT
            old_ref = record.REF
            old_alt = record.ALT[0]

            record.REF = old_alt
            record.ALT = [old_ref]

            # Flip genotype coding
            record.genotypes = [
                [
                    1 if a == 0 else 0 if a == 1 else a,
                    1 if b == 0 else 0 if b == 1 else b,
                    phased
                ]
                for a, b, phased in record.genotypes
            ]

        writer.write_record(record)

    writer.close()
    vcf.close()

    _finalize_output(
        write_path,
        final_path,
        temporary
    )




def ref_to_ancestral(input_vcf, output_vcf):
    """
    Ensure every VCF record has an AA INFO field.

    If AA exists, it is retained.
    If AA is missing, AA is set to REF.

    Parameters
    ----------
    input_vcf : str
        Input VCF(.vcf.gz).
    output_vcf : str
        Output VCF(.vcf.gz).
    """

    write_path, final_path, temporary = _prepare_output(
        input_vcf, output_vcf
    )

    vcf = VCF(input_vcf)

    # Add AA header if absent
    if "ID=AA," not in vcf.raw_header:
        vcf.add_info_to_header({
            "ID": "AA",
            "Description": "Ancestral allele",
            "Type": "String",
            "Number": "1"
        })

    writer = Writer(write_path, vcf)

    for record in vcf:

        aa = record.INFO.get("AA")

        if aa in (None, ".", ""):
            record.INFO["AA"] = record.REF

        writer.write_record(record)

    writer.close()
    vcf.close()

    _finalize_output(
        write_path,
        final_path,
        temporary
    )

