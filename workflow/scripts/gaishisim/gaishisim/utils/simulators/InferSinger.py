from pathlib import Path
import subprocess
from typing import Optional
import os
import pandas as pd
import tskit

from .singer_to_tskit import write_trees
from .tree_functions import compute_eff_popsize_from_diversity
from .vcf_preprocessing import flip_to_ancestral_ref, ref_to_ancestral


class InferSinger:
    """
    Class for applying Singer (including ts-conversion)
    VCF
      |
      | singer_master
      v
    ARG samples
      |
      | convert_to_tskit
      v
    tskit tree sequence (.trees)

    """

    def __init__(
        self,
        vcf_file: str,
        output_prefix: str,
        mutation_rate: float = 1.2e-8,
        start: int = 0,
        end: Optional[int] = None,

        singer_path: str = "resources/singer/releases",

        Ne: Optional[int] = None,

        recombination_map: Optional[str] = None,
        mutation_map: Optional[str] = None,

        ratio: float = 1.0,

        samples = None,
        thin: int = 20,
        polarization: float = 0.99,#0.5,
        n_posterior_samples: int = 1,

        convert_to_trees: bool = True,

        original_ts = None,
        chromosome: str = "1",
        sequence_length = None,
        Ne_estimate_method="diversity", #diversity: from nucleotide diversity, "demes": read from demes, None: default value

        set_flip_to_ancestral_ref = True,

        poplabels_file: Optional[str] = None,
        add_metadata: bool = True,

    ):

        self.vcf_file = Path(vcf_file).resolve()


        # Create final VCF filename if ancestral-reference processing is enabled
        if set_flip_to_ancestral_ref:
            self.vcf_file_final = self.vcf_file.with_name(
                f"{self.vcf_file.stem}_singer{self.vcf_file.suffix}"
            )
        else:
            self.vcf_file_final = self.vcf_file

        # SINGER expects VCF prefix, not filename with .vcf
        if self.vcf_file.suffix == ".vcf":
            self.vcf_prefix = self.vcf_file.with_suffix("")
        else:
            self.vcf_prefix = self.vcf_file

        if self.vcf_file_final.suffix == ".vcf":
            self.vcf_prefix_final = self.vcf_file_final.with_suffix("")
        else:
            self.vcf_prefix_final = self.vcf_file_final


        self.output_prefix = Path(output_prefix).resolve()

        self.output_dir = self.output_prefix.parent
        self.output_name = self.output_prefix.name


        self.mutation_rate = mutation_rate

        self.start = start
        self.end = end

        self.Ne = Ne

        self.recombination_map = recombination_map
        self.mutation_map = mutation_map

        self.ratio = ratio

        #poplabel metadata
        self.poplabels_file = (
            Path(poplabels_file).resolve()
            if poplabels_file is not None
            else None
        )
        self.add_metadata = add_metadata


        self.thin = thin
        self.polarization = polarization


        self.convert_to_trees = convert_to_trees


        self.singer_path = Path(singer_path).resolve()

        self.singer_bin = (
            self.singer_path / "singer_master"
        ).resolve()

        self.convert_bin = (
            self.singer_path / "convert_to_tskit"
        ).resolve()


        self.trees_file = (
            self.output_prefix.with_suffix(".trees")
        )

        self.original_ts = original_ts

        self.chromosome = chromosome

        self.n_posterior_samples = n_posterior_samples

        
        self.samples = samples

        self.sequence_length = sequence_length

        self.Ne_estimate_method = Ne_estimate_method

        self.set_flip_to_ancestral_ref = set_flip_to_ancestral_ref


    # poplabel metadata
    def _add_tree_sequence_metadata(self, ts_file, out_file):
        """
        Add individual and population metadata to a SINGER tree sequence.

        The poplabels file is expected to contain at least:
            sample population group sex 
            population and group can be identica, sex NA

        Samples are assumed to be diploid and therefore represented by
        two nodes per individual.
        """

        if self.poplabels_file is None:
            raise ValueError(
                "poplabels_file must be provided when adding metadata."
            )

        ts_file = Path(ts_file)
        out_file = Path(out_file)

        ts = tskit.load(ts_file)
        tables = ts.dump_tables()

        metadatain = pd.read_csv(
            self.poplabels_file,
            sep=r"\s+"
        )

        required_columns = {"sample", "population"}

        missing = required_columns - set(metadatain.columns)

        if missing:
            raise ValueError(
                f"Missing required columns in poplabels file: {missing}"
            )

        if len(metadatain) * 2 != ts.num_samples:
            raise ValueError(
                f"Number of individuals in poplabels file "
                f"({len(metadatain)}) does not match the number of "
                f"diploid individuals implied by the tree sequence "
                f"({ts.num_samples // 2})."
            )

        tables.individuals.metadata_schema = (
            tskit.MetadataSchema.permissive_json()
        )

        tables.populations.metadata_schema = (
            tskit.MetadataSchema.permissive_json()
        )

        # Create populations
        pop_name_to_id = {}

        for pop_name in metadatain["population"].unique():
            pop_id = tables.populations.add_row(
                metadata={"name": pop_name}
            )
            pop_name_to_id[pop_name] = pop_id

        # update the tables
        # refresh individual tables
        tables.individuals.clear()

        for row in metadatain.itertuples():
            tables.individuals.add_row(
                metadata={"name": row.sample}
            )

        # nodes to indivduals (1 ind, 2 nodes), assumes diploid 
        for ind_id, row in enumerate(metadatain.itertuples()):

            node_1 = 2 * ind_id
            node_2 = node_1 + 1

            target_pop_id = pop_name_to_id[row.population]
            # node 1 of this ind
            tables.nodes[node_1] = tables.nodes[node_1].replace(
                individual=ind_id,
                population=target_pop_id
            )
            # node 1 of this ind
            tables.nodes[node_2] = tables.nodes[node_2].replace(
                individual=ind_id,
                population=target_pop_id
            )
        #output
        ts_wmeta = tables.tree_sequence()

        ts_wmeta.dump(out_file)

        print(
            f"Added metadata to tree sequence: {out_file}"
        )


    # MAIN run function

    def run(self):

        print("Running SINGER")

        if self.set_flip_to_ancestral_ref:

            #if self.ref_to_ancestral:
            ref_to_ancestral(self.vcf_file, self.vcf_file_final)
            flip_to_ancestral_ref(self.vcf_file_final, self.vcf_file_final)

        self._check_binaries()

        self._run_singer()

        if self.convert_to_trees:
            self._convert_to_trees()


    # ======================================================
    # RUN SINGER
    # ======================================================

    def _run_singer(self):

        if self.end is None:
            if self.sequence_length is not None:
                self.end = self.sequence_length
            else:
                self.end = self._get_vcf_end()



        cmd = [

            str(self.singer_bin),

            "-m",
            str(self.mutation_rate),

            "-vcf",
            str(self.vcf_prefix_final),

            "-output",
            str(self.output_prefix),

            "-start",
            str(self.start),

            "-end",
            str(self.end),

            "-n",
            str(self.n_posterior_samples),

            "-thin",
            str(self.thin),

            "-polar",
            str(self.polarization),

        ]



        if self.original_ts is None:
            original_ts = None
        elif isinstance(self.original_ts, (str, Path)):
            original_ts = tskit.load(self.original_ts)
        else:
            original_ts = self.original_ts

        if self.Ne_estimate_method == "diversity":
            self.Ne =  compute_eff_popsize_from_diversity(mutation_rate=self.mutation_rate, vcf_file=self.vcf_file_final, 
                                        ts=original_ts, chr=self.chromosome)


        if self.Ne is not None:
            cmd.extend(
                [
                    "-Ne",
                    str(self.Ne)
                ]
            )


        if self.recombination_map is not None:

            cmd.extend(
                [
                    "-recomb_map",
                    str(self.recombination_map)
                ]
            )


        if self.mutation_map is not None:

            cmd.extend(
                [
                    "-mut_map",
                    str(self.mutation_map)
                ]
            )


        else:

            cmd.extend(
                [
                    "-ratio",
                    str(self.ratio)
                ]
            )


        self._run_command(
            cmd,
            "SINGER inference",
            cwd=self.output_dir
        )


    # ======================================================
    # CONVERT TO tskit
    # ======================================================

    def _get_number_of_singer_samples(self):



        i = 0

        while True:
            node_file = Path(
                f"{self.output_prefix}_nodes_{i}.txt"
            )

            if not node_file.exists():
                break

            i += 1

        if i == 0:
            raise FileNotFoundError(
                f"No SINGER node files found for prefix "
                f"{self.output_prefix}"
            )

        return i

    def _convert_to_trees_default(self):

        if not self.samples:
            arg_samples_to_use = self._get_number_of_singer_samples()
        else:
            arg_samples_to_use = self.samples

    
        cmd = [

            str(self.convert_bin),

            "-input",
            str(self.output_prefix),

            "-output",
            str(self.output_prefix),

            "-start",
            "0",

            "-end",
            str(arg_samples_to_use),

            "-step",
            "1",

        ]

        self._run_command(
            cmd,
            "Convert ARG to tskit",
            cwd=self.output_dir
        )



    def _convert_to_trees(
        self,
        start=0,
        arg_samples_to_use=None,
        step=1,
        overwrite_original_trees_with_metadata=True,
    ):

        if arg_samples_to_use is None:
            if self.samples is None:
                arg_samples_to_use = self._get_number_of_singer_samples()
            else:
                arg_samples_to_use = self.samples

        write_trees(
            input_prefix=self.output_prefix,
            output_prefix=self.output_prefix,
            start=start,
            end=arg_samples_to_use,
            step=step,
        )

        tree_files = []

        for i in range(start, arg_samples_to_use, step):

            tree_file = Path(
                f"{self.output_prefix}_{i}.trees"
            )

            if not tree_file.exists():
                raise FileNotFoundError(
                    f"Expected SINGER tree sequence was not created: "
                    f"{tree_file}"
                )

            if self.add_metadata:

                if self.poplabels_file is None:
                    raise ValueError(
                        "add_metadata=True but no poplabels_file was provided."
                    )

                if overwrite_original_trees_with_metadata:
                    metadata_file = tree_file.with_suffix(
                        ".metadata.tmp.trees"
                    )
                    final_file = tree_file
                else:
                    metadata_file = tree_file.with_name(
                        f"{tree_file.stem}.metadata.trees"
                    )
                    final_file = metadata_file

                self._add_tree_sequence_metadata(
                    ts_file=tree_file,
                    out_file=metadata_file,
                )

                if overwrite_original_trees_with_metadata:
                    os.replace(metadata_file, final_file)

            else:
                final_file = tree_file

            tree_files.append(final_file)

        return tree_files



    #check existing singer binaries
    def _check_binaries(self):

        if not self.singer_bin.exists():
            raise FileNotFoundError(
                f"SINGER executable missing: {self.singer_bin}"
            )


        if self.convert_to_trees:

            if not self.convert_bin.exists():
                raise FileNotFoundError(
                    f"SINGER converter missing: {self.convert_bin}"
                )
            
    def _get_vcf_end(self):
        """
        Return the genomic position of the last SNP in the VCF.
        Assumes a standard VCF file.
        """
        last_pos = None

        with open(self.vcf_file_final, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.split("\t")
                last_pos = int(fields[1])  # POS column

        if last_pos is None:
            raise ValueError("No variants found in VCF")

        return last_pos


    # ======================================================
    # EXECUTE
    # ======================================================

    def _run_command(
        self,
        cmd,
        label,
        cwd=None
    ):

        print(f"\n[{label}]")
        print(" ".join(map(str, cmd)))


        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            cwd=cwd,
        )


        if result.returncode != 0:

            print(result.stderr)

            raise RuntimeError(
                f"{label} failed"
            )


        print("InferSinger process finished!")