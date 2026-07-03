import os
import subprocess
from pathlib import Path
from typing import Optional
import demes

from .simulation_utils import effective_population_size


class InferRelate:
    """
    Full Relate pipeline with:

    msprime VCF → Relate inference → RelateParallel (optional)
    → tskit tree sequence conversion

    Ne policy:
        1. demes (true simulation Ne)
        2. fixed Ne
        3. Relate estimation
    """

    def __init__(
        self,
        vcf_file: str,
        output_prefix: str,
        recombination_rate: float,
        mutation_rate: float,
        sequence_length: int = None,
        relate_path: str = "resources/relate",
        genetic_map: Optional[str] = None,
        demes_file: Optional[str] = None,
        Ne: Optional[int] = None,
        mode: str = "standard",  # standard | parallel
        chromosome: str = "1",
        threads: int = 1,
        memory: int = 5,
        try_to_estimate_Ne_from_demes=False
    ):
        self.vcf_file = Path(vcf_file)
        #self.output_prefix = Path(output_prefix)
        self.output_prefix = Path(output_prefix).resolve()

        # directory where outputs should go
        self.output_dir = self.output_prefix.parent

        # name for RelateParallel 
        self.relate_output_name = self.output_prefix.name

        self.recombination_rate = recombination_rate
        self.sequence_length = sequence_length
        self.mutation_rate = mutation_rate

        #self.relate_path = Path(relate_path)
        self.genetic_map = genetic_map
        self.demes_file = demes_file
        self.Ne = Ne

        self.mode = mode
        self.chromosome = chromosome
        self.threads = threads
        self.memory = memory
        self.try_to_estimate_Ne_from_demes = try_to_estimate_Ne_from_demes


        self.relate_path = Path(relate_path).resolve()

        self.relate_bin = (self.relate_path / "bin" / "Relate").resolve()
        self.convert_bin = (self.relate_path / "bin" / "RelateFileFormats").resolve()

        self.parallel_script = (
            self.relate_path / "scripts" / "RelateParallel" / "RelateParallel.sh"
        ).resolve()

        '''
        # binaries
        self.relate_bin = self.relate_path / "bin" / "Relate"
        self.convert_bin = self.relate_path / "bin" / "RelateFileFormats"
        self.parallel_script = (
            self.relate_path / "scripts" / "RelateParallel" / "RelateParallel.sh"
        )
        

        # files
        self.haps_file = f"{self.output_prefix}.haps"
        self.sample_file = f"{self.output_prefix}.sample"
        self.map_file = (
            self.genetic_map if self.genetic_map else f"{self.output_prefix}.map"
        )

        
        self.output_prefix = Path(output_prefix).resolve()
        self.output_dir = self.output_prefix.parent
        self.relate_output_name = self.output_prefix.name
        '''
        self.haps_file = (self.output_dir / f"{self.relate_output_name}.haps").resolve()
        self.sample_file = (self.output_dir / f"{self.relate_output_name}.sample").resolve()
        self.map_file = (self.output_dir / f"{self.relate_output_name}.map").resolve()


        self.ne_script = self.relate_path / "bin" / "EstimatePopulationSize.sh"
        self.popsize_file = f"{self.output_prefix}.PopSize.estimated_pop_size.txt"

        self.trees_file = f"{self.output_prefix}.trees"

    # ==========================================================
    # MAIN PIPELINE
    # ==========================================================

    def run(self):
        print("InferRelate new version")

        self._check_binaries()

        self._convert_vcf()
        self._ensure_map()

        Ne_value = self._resolve_Ne()

        if self.mode == "parallel":
            self._run_relate_parallel(Ne_value)
        else:
            self._run_relate_standard(Ne_value)

        self._convert_to_trees()

    # ==========================================================
    # NE POLICY
    # ==========================================================


    def _resolve_Ne(self) -> int:
        if self.Ne is not None:
            return self.Ne

        if self.demes_file is not None and self.try_to_estimate_Ne_from_demes:
            return self._Ne_from_demes()

        #default value
        return 30000
        
    

    def _Ne_from_demes(self) -> int:
        """
        Compute a single Ne value from demes YAML file.
        Uses recent extreme by default.
        """

        ne_dict = effective_population_size(
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
            "-N", str(Ne),

            # IMPORTANT: basename only
            "-o", self.relate_output_name,

            "--threads", str(self.threads),
        ]

        # Run in correct output directory so files are written properly
        self._run_command(cmd, "RelateParallel", cwd=self.output_dir)

    # ==========================================================
    # CONVERT TO tskit TREE SEQUENCE
    # ==========================================================



    def _convert_to_trees(self):

        prefix = str(self.output_prefix)

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
        total_bp = self.sequence_length
        total_cM = self.sequence_length * self.recombination_rate * 100.0

        with open(self.map_file, "w") as f:
            # position, recombination rate (cM/Mb), genetic position (cM)
            
            f.write(f"0\t{self.recombination_rate * 100.0}\t0.0\n")
            
            f.write(f"{total_bp}\t{self.recombination_rate * 100.0}\t{total_cM}\n")

    # ==========================================================
    # workflow functions
    # ==========================================================

    def _check_binaries(self):
        if not self.relate_bin.exists():
            raise FileNotFoundError("Relate binary missing")
        if not self.convert_bin.exists():
            raise FileNotFoundError("RelateFileFormats missing")
        if self.mode == "parallel" and not self.parallel_script.exists():
            raise FileNotFoundError("RelateParallel missing")


        
    def _run_command(self, cmd, label, cwd=None):
        print(f"\n[{label}]")
        print(" ".join(map(str, cmd)))

        result = subprocess.run(
            cmd,
            shell=False,
            capture_output=True,
            text=True,
            cwd=cwd  
        )

        if result.returncode != 0:
            print(result.stderr)
            raise RuntimeError(f"{label} failed")

        print("done")
