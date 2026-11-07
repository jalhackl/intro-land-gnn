import os
import subprocess
from pathlib import Path
from typing import Optional
import demes
import tskit

import json
import tszip

from .simulation_utils import effective_population_size_from_yaml
from .tree_functions import compute_eff_popsize_from_diversity


class InferRelate:
    """
    Full Relate pipeline with:

    msprime VCF → Relate inference → RelateParallel (optional)
    → tskit tree sequence conversion

    for calculating Ne, per default diversity is calculated from the vcf-file
    """

    def __init__(
        self,
        vcf_file: str,
        output_prefix: str,
        recombination_rate: float,
        mutation_rate: float =1.2e-8,
        sequence_length: int = None,
        relate_path: str = "resources/relate",
        genetic_map: Optional[str] = None,
        demes_file: Optional[str] = None,
        Ne: Optional[int] = None,
        mode: str = "standard",  # standard | parallel
        chromosome: str = "1",
        threads: int = 1,
        memory: int = 5,
        Ne_estimate_method="diversity", #diversity: from nucleotide diversity, "demes": read from demes, None: default value (usually diversity is recommended)
        default_Ne: int = 30000,
        original_ts_file = None,
        create_relate_poplabels: bool = True,
        # alternatively, it should be possible to load a poplabels file
        #poplabels_file = None,
        add_poplabels_metadata: bool = True,
        reestimate_branch_lengths: bool = True,
        overwrite_original_trees_with_reestimated=False,
        also_process_not_reestimated_versions=True,
        years_per_gen=29
    ):
        self.years_per_gen = years_per_gen
        
        self.vcf_file = Path(vcf_file)
        self.output_prefix = Path(output_prefix).resolve()

        # directory where outputs should go
        self.output_dir = self.output_prefix.parent

        # name for RelateParallel 
        self.relate_output_name = self.output_prefix.name


        # name for branch length updated version
        # (application of reestimatepopulationsizes)
        self.output_prefix_reestimated = (
            self.output_prefix.parent
            / (self.output_prefix.name + "_reestimated")
        ).resolve()
        self.relate_output_name_reestimated = self.output_prefix.name + "_reestimated"

        self.recombination_rate = float(recombination_rate)
        self.sequence_length = (
            int(sequence_length) if sequence_length is not None else None
        )
        self.mutation_rate = float(mutation_rate)

        #self.relate_path = Path(relate_path)
        self.genetic_map = genetic_map
        self.demes_file = demes_file
        self.Ne = Ne

        self.mode = mode
        self.chromosome = chromosome
        self.threads = threads
        self.memory = memory
        self.Ne_estimate_method = Ne_estimate_method


        self.relate_path = Path(relate_path).resolve()

        self.relate_bin = (self.relate_path / "bin" / "Relate").resolve()
        self.convert_bin = (self.relate_path / "bin" / "RelateFileFormats").resolve()

        self.parallel_script = (
            self.relate_path / "scripts" / "RelateParallel" / "RelateParallel.sh"
        ).resolve()


        self.haps_file = (self.output_dir / f"{self.relate_output_name}.haps").resolve()
        self.sample_file = (self.output_dir / f"{self.relate_output_name}.sample").resolve()
        self.map_file = (self.output_dir / f"{self.relate_output_name}.map").resolve()


        self.ne_script = self.relate_path / "scripts" / "EstimatePopulationSize" / "EstimatePopulationSize.sh"
        self.reestimate_branch_lengths = reestimate_branch_lengths
        self.overwrite_original_trees_with_reestimated = overwrite_original_trees_with_reestimated
        self.also_process_not_reestimated_versions = also_process_not_reestimated_versions


        self.popsize_file = f"{self.output_prefix}.PopSize.estimated_pop_size.txt"

        self.trees_file = f"{self.output_prefix}.trees"
        #reestimated version
        self.trees_file_reestimated = f"{self.output_prefix_reestimated}.trees"
        
        
        self.default_Ne = default_Ne
        self.original_ts_file = original_ts_file

        #poplabel file - necessary for downstream analyses
        self.create_relate_poplabels = create_relate_poplabels
        self.add_poplabels_metadata = add_poplabels_metadata

        self.poplabels_file = (
            self.output_dir / f"{self.relate_output_name}.relate.poplabels"
        ).resolve()
        #duplicate!
        #self.poplabels_output_file = self.output_dir / f"{self.relate_output_name}.relate.poplabels"


        self.metadata_trees_file = (
            self.output_dir / f"{self.relate_output_name}.metadata.trees"
        ).resolve()
        #reestimated version
        self.metadata_trees_file_reestimated = (
            self.output_dir / f"{self.relate_output_name_reestimated}.metadata.trees"
        ).resolve()

    # ==========================================================
    # MAIN PIPELINE
    # ==========================================================


    def run(self):

        self._check_binaries()

        self._convert_vcf()
        self._ensure_map()

        Ne_value = self._resolve_Ne()

        if self.mode == "parallel":
            self._run_relate_parallel(Ne_value)
        else:
            self._run_relate_standard(Ne_value)

        if self.create_relate_poplabels:
            self._create_relate_poplabels()

        # call estimatepopulationsizes script
        if self.reestimate_branch_lengths: 
            self._run_ne_Script()

            # just the trees file?
            #if self.overwrite_original_trees_with_reestimated:
            #    os.replace(self.relate_, final_file)

            self._convert_to_trees(self.output_prefix_reestimated)

            if self.add_poplabels_metadata:
                self._add_poplabels_metadata(ts_file=self.trees_file_reestimated,
                                             out_trees_file=self.metadata_trees_file_reestimated)


        #overwrite_original_trees_with_reestimated = self.overwrite_original_trees_with_reestimated

        if self.also_process_not_reestimated_versions:

            self._convert_to_trees()

            if self.add_poplabels_metadata:
                self._add_poplabels_metadata()



    def _run_ne_Script(self):
        cmd = [
            str(self.ne_script),
            "--input",  self.relate_output_name,
            "--poplabels", self.poplabels_file,
            "-m", str(self.mutation_rate),
            "--years_per_gen", str(self.years_per_gen), 
            "-o", self.relate_output_name_reestimated,
            "--noplot" 
        ]

        self._run_command(cmd, "Relate reestimatiaton Popsizes script", cwd=self.output_dir)    

    ## optional arguments according to relate documentation
    '''
    --noanc
        Optional: Flag specifying whether to reinfer branch lengths at the end. 0: reinfer, 1: don't reinfer. Default: 0.
    --pop_of_interest
        Optional: Specifies populations for which population size is estimated. The second column of the .poplabels file is used for population labels. You can specify multiple populations separated by commas, e.g. pop1,pop2. When using this option, it is recommended (but NOT necessary), to run Relate with --annot. This can also be done after running Relate without --annot, see here. 
    --first_chr & --last_chr
        Optional: Allow for multiple input files, where filenames are assumed to be numbered by chromosome index.
        Example: -i example --first_chr 1 --last_chr 2 assumes example_chr1.anc, example_chr1.mut, example_chr2.anc, example_chr2.mut as input.
    --chr
        Optional: File containing chromosome names (one per line). Overrides --first_chr and --last_chr.
        Example: -i example --chr chr.txt assumes example_chr*.anc, example_chr*.mut as input.
    --threads
        Optional: Maximum number of threads. Default is 1.
    --bins
        Optional: Specify epoch bins. Format: lower, upper, stepsize for function c(0,10^seq(lower, upper, stepsize)). 
    --years_per_gen
        Optional: Years per generation. Default: 28.0
    --num_iter
        Optional: Number of iterations of the algorithm. Default is 5.
    --threshold
        Optional: Float in [0,1] which is the fraction of trees to be dropped, ranked by number of mutations mapping. This option serves two purposes: First, it removes bad/uninformative trees. Second it speeds up the inference. As we increase this value, we remove more trees. Default: 0.5.
    --seed
        Optional: Seed for the random number generator used in the MCMC algorithm for estimating branch lengths.
    --noplot
        Optional: Takes no value. If specified, do not plot the results using R.
    '''


    # ==========================================================
    # NE POLICY
    # ==========================================================


    def _resolve_Ne(self) -> int:
        if self.Ne is not None:
            return self.Ne
        
        if self.Ne_estimate_method == "diversity":

            from pathlib import Path
            import tskit

            if self.original_ts_file is None:
                original_ts = None
            elif isinstance(self.original_ts_file, (str, Path)):
                original_ts = tskit.load(self.original_ts_file)
            else:
                original_ts = self.original_ts_file

            return compute_eff_popsize_from_diversity(
                mutation_rate=self.mutation_rate,
                vcf_file=self.vcf_file,
                ts=original_ts,
                chr=self.chromosome
            )



        if self.demes_file is not None and self.Ne_estimate_method == "demes":
            return self._Ne_from_demes()

        #default value
        return self.default_Ne
        
    

    def _Ne_from_demes(self) -> int:
        """
        Compute a single Ne value from demes YAML file.
        Uses recent extreme by default.
        """

        ne_dict = effective_population_size_from_yaml(
            yaml_path=self.demes_file,
            mode="recent",
            summary="extreme",
        )

        if not ne_dict:
            raise ValueError("No Ne values computed from demes file")

        # If multiple demes returned, average them
        values = [v for v in ne_dict.values() if v is not None]

        if not values:
            raise ValueError("Computed Ne values are empty")

        Ne_value = int(sum(values) / len(values))

        print(f"[Ne from demes] Using Ne = {Ne_value}")

        return Ne_value

        

    # ==========================================================
    # RELATE STANDARD
    # ==========================================================

    def _run_relate_standard(self, Ne: int):
        cmd = [
            str(self.relate_bin),
            "--mode", "All",
            "--haps", self.haps_file,
            "--sample", self.sample_file,
            "--map", self.map_file,
            "-N", str(Ne),
            "-m", str(self.mutation_rate),
            "--memory", str(self.memory),
            "-o", self.relate_output_name,  
        ]

        self._run_command(cmd, "Relate standard", cwd=self.output_dir)    

    # ==========================================================
    # RELATE PARALLEL
    # ==========================================================

    def _run_relate_parallel(self, Ne: int):
        print("the haps file2")
        print(self.haps_file)
        print(self.parallel_script)

        cmd = [
            #"bash",
            str(self.parallel_script),

            "--haps", self.haps_file,
            "--sample", self.sample_file,
            "--map", self.map_file,

            "-m", str(self.mutation_rate),
            "--memory", str(self.memory),
            "-N", str(Ne),

            # basename only
            "-o", self.relate_output_name,

            "--threads", str(self.threads),
        ]

        # Run in correct output directory so files are written properly
        self._run_command(cmd, "RelateParallel", cwd=self.output_dir)

    # CONVERT TO tskit tree sequence command

    def _convert_to_trees(self, prefix=None):

        prefix = str(prefix if prefix is not None else self.output_prefix)

        cmd = [
            str(self.convert_bin),
            "--mode", "ConvertToTreeSequence",
            "-i", prefix,
            "-o", prefix,
        ]

        self._run_command(cmd, "Convert to tskit trees")

    # ==========================================================
    # VCF CONVERSION
    # ==========================================================

    def _convert_vcf(self):
        cmd = [
            str(self.convert_bin),
            "--mode",
            "ConvertFromVcf",
            "--input",
            str(self.vcf_file.with_suffix("")),
            "--haps",
            self.haps_file,
            "--sample",
            self.sample_file,
        ]
        self._run_command(cmd, "VCF conversion")

    # ==========================================================
    # MAP
    # ==========================================================

    def _ensure_map(self):
        if self.genetic_map:
            return
        self._create_uniform_map()



    def _create_uniform_map(self):
    
        if self.sequence_length is None:
            raise ValueError(
                "sequence_length must be provided when no genetic map is supplied."
            )
    
        total_bp = self.sequence_length
        total_cM = self.sequence_length * self.recombination_rate * 100.0

        with open(self.map_file, "w") as f:
            # position, recombination rate (cM/Mb), genetic position (cM)
            
            f.write(f"0\t{self.recombination_rate * 100.0}\t0.0\n")
            
            f.write(f"{total_bp}\t{self.recombination_rate * 100.0}\t{total_cM}\n")

    # ==========================================================
    # workflow functions
    # ==========================================================

    # check existence of relate binaries
    def _check_binaries(self):
        if not self.relate_bin.exists():
            raise FileNotFoundError("Relate binary missing")
        if not self.convert_bin.exists():
            raise FileNotFoundError("RelateFileFormats missing")
        if self.mode == "parallel" and not self.parallel_script.exists():
            raise FileNotFoundError("RelateParallel missing")


        

    def _run_command(self, cmd, label, cwd=None):
        print(f"\n[{label}]")
        print(f"Working directory: {cwd}")
        print("Command:", " ".join(map(str, cmd)), flush=True)

        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            cwd=cwd,
        )

        if result.stdout:
            print(f"\n[{label} stdout]")
            print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")

        if result.stderr:
            print(f"\n[{label} stderr]")
            print(result.stderr, end="" if result.stderr.endswith("\n") else "\n")

        print(f"\n[{label}] Return code: {result.returncode}", flush=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"{label} failed with exit code {result.returncode}"
            )

        print(f"{label} finished successfully!", flush=True)



    def _write_relate_poplabels(
        self,
        ts: tskit.TreeSequence,
        sample_file: str,
        output_file: str,
    ) -> None:
        """
        Write a Relate-compatible .poplabels file. (sample population group sex)
        The samples are written in the same order as the Relate
        .sample file.

        msprime sample IDs are expected to be tsk_<individual_id>
        for other input, a different function has to be used

        Format:
            sample population group sex
        """

        with open(sample_file) as f:
            lines = [line.strip().split() for line in f if line.strip()]

        # Relate .sample has a header on the first line.
        # and also the second line has to be omitted
        samples = lines[2:]

        with open(output_file, "w") as out:
            out.write("sample population group sex\n")

            for row in samples:
                id1 = row[0]
                id2 = row[1]

                if id2 == "NA":
                    # Haploid
                    sample_id = id1
                    sex = 1
                else:
                    # Diploid
                    sample_id = id1
                    sex = "NA"

                # msprime sample IDs are expected to be tsk_<individual_id>
                # for other input, a different function has to be used
                try:
                    individual_id = int(sample_id.split("_")[-1])
                except ValueError as e:
                    raise ValueError(
                        f"Could not extract individual ID from '{sample_id}'. "
                        "Expected an ID such as 'tsk_0'."
                    ) from e

                if individual_id >= ts.num_individuals:
                    raise ValueError(
                        f"Individual {individual_id} from '{sample_id}' "
                        f"does not exist in the tree sequence "
                        f"(num_individuals={ts.num_individuals})."
                    )

                individual = ts.individual(individual_id)

                if len(individual.nodes) == 0:
                    raise ValueError(
                        f"Individual {individual_id} ('{sample_id}') "
                        "has no nodes in the tree sequence."
                    )

                # use first node to determine the population.
                node = ts.node(individual.nodes[0])
                population = ts.population(node.population)

                pop_label = population.metadata.get(
                    "name",
                    str(node.population),
                )

                group_label = pop_label

                out.write(
                    f"{sample_id} {pop_label} {group_label} {sex}\n"
                )



    def _create_relate_poplabels(self):
        if self.original_ts_file is None:
            raise ValueError(
                "create_relate_poplabels=True currently requires "
                "original_ts_file to be provided (of course, one could perhaps also load the data from the ind-files and/or the vcf)."
            )


        ts = tskit.load(self.original_ts_file)


        #output_file = self.output_dir / f"{self.relate_output_name}.relate.poplabels"

        self._write_relate_poplabels(
            ts=ts,
            sample_file=str(self.sample_file),
            output_file=str(self.poplabels_file),
        )

        print(f"Created Relate poplabels: {self.poplabels_file}")


    def _add_poplabels_metadata_tszip(self):

        if not self.poplabels_file.exists():
            raise FileNotFoundError(
                f"Poplabels file not found: {self.poplabels_file}"
            )

        if not Path(self.trees_file).exists():
            raise FileNotFoundError(
                f"Relate tree sequence not found: {self.trees_file}"
            )

        print(f"Loading Relate tree sequence from {self.trees_file}")
        ts = tszip.load(self.trees_file)

        tables = ts.dump_tables()

        # --------------------------------------------------------------
        # Read poplabels
        # --------------------------------------------------------------

        print(
            f"Loading population and individual information from "
            f"{self.poplabels_file}"
        )

        with open(self.poplabels_file) as f:
            header = f.readline().strip().split()

            expected_header = ["sample", "population", "group", "sex"]

            if header != expected_header:
                raise ValueError(
                    f"Unexpected poplabels header: {header}. "
                    f"Expected: {expected_header}"
                )

            metadata_rows = []

            for line_number, line in enumerate(f, start=2):
                fields = line.strip().split()

                if len(fields) != len(header):
                    raise ValueError(
                        f"Malformed poplabels line {line_number}: {fields}"
                    )

                metadata_rows.append(dict(zip(header, fields)))

        print(f"- Number of individuals: {len(metadata_rows)}")

        # --------------------------------------------------------------
        # Populations
        # --------------------------------------------------------------

        populations = sorted(
            {row["population"] for row in metadata_rows}
        )

        print(f"- Number of unique populations: {len(populations)}")

        population_to_id = {
            pop: i for i, pop in enumerate(populations)
        }

        # --------------------------------------------------------------
        # Validate individual count
        # --------------------------------------------------------------

        if len(metadata_rows) != tables.individuals.num_rows:
            raise ValueError(
                "Number of individuals in .poplabels does not match "
                f"the Relate tree sequence: "
                f"{len(metadata_rows)} poplabels vs "
                f"{tables.individuals.num_rows} individuals."
            )

        # --------------------------------------------------------------
        # Individual metadata
        # --------------------------------------------------------------

        print("Setting individual metadata")

        tables.individuals.metadata_schema = tskit.MetadataSchema(None)

        metadata_column = [
            json.dumps(row).encode("utf-8")
            for row in metadata_rows
        ]

        md, md_offset = tskit.pack_bytes(metadata_column)

        tables.individuals.set_columns(
            flags=tables.individuals.flags,
            location=tables.individuals.location,
            location_offset=tables.individuals.location_offset,
            parents=tables.individuals.parents,
            parents_offset=tables.individuals.parents_offset,
            metadata=md,
            metadata_offset=md_offset,
        )

        # --------------------------------------------------------------
        # Population metadata
        # --------------------------------------------------------------

        print("Setting population metadata")

        tables.populations.metadata_schema = tskit.MetadataSchema(None)

        population_metadata = [
            json.dumps({"name": pop}).encode("utf-8")
            for pop in populations
        ]

        md, md_offset = tskit.pack_bytes(population_metadata)

        tables.populations.set_columns(
            metadata=md,
            metadata_offset=md_offset,
        )

        # --------------------------------------------------------------
        # Set individual IDs on sample nodes
        # --------------------------------------------------------------

        print("Setting individual IDs on sample nodes")

        sample_nodes = ts.samples()

        if len(sample_nodes) != 2 * len(metadata_rows):
            raise ValueError(
                "Unexpected number of sample nodes: "
                f"{len(sample_nodes)} nodes for "
                f"{len(metadata_rows)} individuals. "
                "This implementation expects diploid Relate samples."
            )

        nodes_individual = tables.nodes.individual.copy()

        for individual_id, node_id_pair in enumerate(
            sample_nodes.reshape((-1, 2))
        ):
            for node_id in node_id_pair:
                nodes_individual[node_id] = individual_id

        tables.nodes.individual = nodes_individual

        # --------------------------------------------------------------
        # Set population IDs on sample nodes
        # --------------------------------------------------------------

        print("Setting population IDs on sample nodes")

        nodes_population = tables.nodes.population.copy()

        for individual_id, row in enumerate(metadata_rows):

            population_id = population_to_id[row["population"]]

            for node_id in ts.individual(individual_id).nodes:
                nodes_population[node_id] = population_id

        tables.nodes.population = nodes_population

        # Create new tree sequence
        new_ts = tables.tree_sequence()

        tszip.compress(
            new_ts,
            str(self.metadata_trees_file),
        )

        print(
            f"Done. Metadata-enriched tree sequence: "
            f"{self.metadata_trees_file}"
        )


    def _add_poplabels_metadata(self, ts_file=None, poplabels_file=None, out_trees_file=None):
        ts_file = str(ts_file if ts_file is not None else self.trees_file)
        poplabels_file = str(poplabels_file if poplabels_file is not None else self.poplabels_file)
        out_trees_file = str(out_trees_file if out_trees_file is not None else self.metadata_trees_file)


        print(f"ts_file: {ts_file}")
        print(f"poplabels_file: {poplabels_file}")
        print(f"out_trees_file: {out_trees_file}")

        print(f"Loading tree sequence from {ts_file}")
        ts = tskit.load(ts_file)

        tables = ts.dump_tables()

        # load .poplabels
        print(f"Loading population and individual info from {poplabels_file}")

        uniq_pops = []
        metadata_column = []

        with open(poplabels_file, "r") as ifl:
            header = ifl.readline().strip().split()

            assert header == ["sample", "population", "group", "sex"]

            for l in ifl:
                line = l.strip().split()

                metadata = {}

                for k, v in zip(header, line):
                    metadata[k] = v

                metadata_column.append(
                    json.dumps(metadata).encode("utf-8")
                )

                uniq_pops.append(line[1])

        uniq_pops = sorted(set(uniq_pops))

        print(f"- Number of individuals: {len(metadata_column)}")
        print(f"- Number of unique populations: {len(uniq_pops)}")

        # set individuals table
        print("Setting individuals table")

        basic_schema = tskit.MetadataSchema(None)
        tables.individuals.metadata_schema = basic_schema

        flags = [0 for x in range(len(metadata_column))]

        encoded_metadata_column = [
            tables.individuals.metadata_schema.validate_and_encode_row(r)
            for r in metadata_column
        ]

        md, md_offset = tskit.pack_bytes(encoded_metadata_column)

        tables.individuals.set_columns(
            flags=flags,
            metadata=md,
            metadata_offset=md_offset,
        )

        # set populations table
        print("Setting populations table")

        population_metadata_column = []

        for pop in uniq_pops:
            population_metadata_column.append(
                json.dumps({"name": pop}).encode("utf-8")
            )

        tables.populations.metadata_schema = basic_schema

        encoded_pop_metadata_column = [
            tables.populations.metadata_schema.validate_and_encode_row(r)
            for r in population_metadata_column
        ]

        md, md_offset = tskit.pack_bytes(encoded_pop_metadata_column)

        tables.populations.set_columns(
            metadata=md,
            metadata_offset=md_offset,
        )

        # set individual id of sample nodes in node table
        print("Setting individual id of sample nodes in node table")

        num_individuals = len(metadata_column)

        new_individual_of_sample_nodes = []

        for i in range(num_individuals):
            new_individual_of_sample_nodes.extend([i, i])

        new_individual = tables.nodes.individual

        new_individual[:len(new_individual_of_sample_nodes)] = (
            new_individual_of_sample_nodes
        )

        tables.nodes.individual = new_individual

        # set population of sample nodes in node table
        print("Setting population id of sample nodes in node table")

        num_indivs = len(tables.individuals)

        new_population_for_sample_nodes = []

        for idx in range(num_indivs):
            indiv = tables.individuals[idx]

            pop = json.loads(
                indiv.metadata.decode("utf-8")
            )["population"]

            pop_idx = uniq_pops.index(pop)

            new_population_for_sample_nodes.extend(
                [pop_idx, pop_idx]
            )

        new_population = tables.nodes.population

        new_population[:len(new_population_for_sample_nodes)] = (
            new_population_for_sample_nodes
        )

        tables.nodes.population = new_population

        # convert to new tree sequence
        new_ts = tables.tree_sequence()

        new_ts.dump(out_trees_file)

        print(f"Done. out tree sequence: {out_trees_file}")