import numpy as np
import pandas as pd


def create_map_file(x, y, map_file="sim.map"):
    with open(map_file, "w") as f:
        f.write("1\t.\t0.0\t0\n")
        f.write(f"1\t.\t{x}\t{y}\n")


def process_sprime_output(in_file, out_file):
    """
    Description:
        Helper function for converting output from SPrime to BED format.

    Arguments:
        in_file str: Name of the input file.
        out_file str: Name of the output file.
    """
    # read in the data - the SPrime output
    df = pd.read_csv(in_file, delimiter="\t")

    # drop columns that are not needed
    df2 = df.drop(["ID", "REF", "ALT", "ALLELE"], axis=1)

    # add columns START and END with the highest and lowest position of each chromosome, segment and score
    df2["START"] = df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].transform(min)
    df2["END"] = df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].transform(max)

    # group by chromosome, segment and score - drop the column position
    df3 = df2.loc[df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].idxmax()]
    df4 = df3.drop(["POS"], axis=1)

    # get the right order (for the bed file)
    df_final = df4[["CHROM", "START", "END", "SEGMENT", "SCORE"]].sort_values(
        by=["START", "SEGMENT"]
    )

    np.savetxt(out_file, df_final.values, fmt="%s", delimiter="\t")

# -----------------------------------------------------------------------------------------------------------------------
# process_sprime_output fn from https://github.com/admixVIE/sstar-analysis/blob/main/utils/utils.py
# -----------------------------------------------------------------------------------------------------------------------

# process_sprime_output fn moved to workflow/scripts/sstar_sprime_out_to_bed.py
#def process_sprime_output(in_file, out_file):
#    """
#    Description:
#        Helper function for converting output from SPrime to BED format.
#
#    Arguments:
#        in_file str: Name of the input file.
#        out_file str: Name of the output file.
#    """
#    # read in the data - the SPrime output
#    df = pd.read_csv(in_file, delimiter="\t")
#
#    # drop columns that are not needed
#    df2 = df.drop(["ID", "REF", "ALT", "ALLELE"], axis=1)
#
#    # add columns START and END with the highest ans lowest position of each chromosome, segment and score
#    df2["START"] = df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].transform(min)
#    df2["END"] = df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].transform(max)
#
#    # group by chromosome, segment and score - drop the column position
#    df3 = df2.loc[df2.groupby(["CHROM", "SCORE", "SEGMENT"])["POS"].idxmax()]
#    df4 = df3.drop(["POS"], axis=1)
#
#    # get the right order (for the bed file)
#    df_final = df4[["CHROM", "START", "END", "SEGMENT", "SCORE"]].sort_values(
#        by=["START", "SEGMENT"]
#    )
#
#    np.savetxt(out_file, df_final.values, fmt="%s", delimiter="\t")
