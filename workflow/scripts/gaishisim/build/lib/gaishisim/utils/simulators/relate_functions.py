import os
import subprocess
import tskit

relate_dir = os.path.join("resources", "relate")
convert_script = os.path.join(relate_dir, "scripts", "ConvertFromVcf.sh")

print(relate_dir)



def _create_relate_files(
    #self,
    ts: tskit.TreeSequence,
    vcf_file: str,
    output_prefix: str,
    output_dir: str,
    relate_dir: str,
) -> None:
    """
    Runs full Relate input workflow:
    1) Convert VCF -> .haps/.sample
    2) Create .poplabels consistent with .sample
    """

    convert_script = os.path.join(relate_dir, "scripts", "ConvertFromVcf.sh")

    relate_prefix = os.path.join(output_dir, output_prefix)

    # ---- 1) Convert VCF to haps/sample ----
    cmd = [
        convert_script,
        "--mode", "All",
        "--vcf", vcf_file,
        "--haps", relate_prefix + ".haps",
        "--sample", relate_prefix + ".sample",
    ]

    subprocess.run(cmd, check=True)

    # ---- 2) Create poplabels ----
    sample_file = relate_prefix + ".sample"
    poplabels_file = relate_prefix + ".poplabels"

    #self.
    _write_relate_poplabels(
        ts=ts,
        sample_file=sample_file,
        output_file=poplabels_file,
    )


def _write_relate_poplabels(
    #self,
    ts: tskit.TreeSequence,
    sample_file: str,
    output_file: str,
) -> None:
    """
    Write Relate-compliant .poplabels file.

    Format:
    sample population group sex
    """

    with open(sample_file) as f:
        lines = [l.strip().split() for l in f]

    header = lines[0]
    samples = lines[1:]  # skip header

    with open(output_file, "w") as out:
        out.write("sample population group sex\n")

        for row in samples:
            id1 = row[0]
            id2 = row[1]

            if id2 == "NA":
                # Haploid case
                sample_id = id1
                sex = 1
            else:
                # Diploid case
                sample_id = id1
                sex = "NA"

            # Extract node index from msprime naming: tsk_{node}
            node_idx = int(sample_id.split("_")[-1])

            node = ts.node(node_idx)
            pop = ts.population(node.population)
            pop_label = pop.metadata.get("name", str(node.population))

            group_label = pop_label

            out.write(f"{sample_id} {pop_label} {group_label} {sex}\n")