import os
import subprocess
from pathlib import Path
from typing import Optional

import demes


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
        sequence_length: int,
        mutation_rate: float,
        relate_path: str = "resources/relate",
        genetic_map: Optional[str] = None,
        demes_file: Optional[str] = None,
        Ne: Optional[int] = None,
        mode: str = "standard",  # standard | parallel
        chromosome: str = "1",
        threads: int = 1,
        memory: int = 5,
    ):
        self.vcf_file = Path(vcf_file)
        self.output_prefix = Path(output_prefix)

        self.recombination_rate = recombination_rate
        self.sequence_length = sequence_length
        self.mutation_rate = mutation_rate

        self.relate_path = Path(relate_path)
        self.genetic_map = genetic_map
        self.demes_file = demes_file
        self.Ne = Ne

        self.mode = mode
        self.chromosome = chromosome
        self.threads = threads
        self.memory = memory

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

        self.trees_file = f"{self.output_prefix}.trees"

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

        self._convert_to_trees()

    # ==========================================================
    # NE POLICY
    # ==========================================================

    def _resolve_Ne(self) -> int:
        if self.demes_file is not None:
            return self._Ne_from_demes()

        if self.Ne is not None:
            return self.Ne

        return self._estimate_Ne()

    def _Ne_from_demes(self) -> int:
        graph = demes.load(self.demes_file)

        sizes = []
        for pop in graph.demes.values():
            for epoch in pop.epochs:
                if epoch.start_size is not None:
                    sizes.append(epoch.start_size)

        if not sizes:
            raise ValueError("No Ne found in demes file")

        return int(sum(sizes) / len(sizes))

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
            "--threads", str(self.threads),
            "--memory", str(self.memory),
            "-o", str(self.output_prefix),
        ]
        self._run_command(cmd, "Relate standard")

    # ==========================================================
    # RELATE PARALLEL
    # ==========================================================

    def _run_relate_parallel(self, Ne: int):
        cmd = [
            str(self.parallel_script),
            "-i", self.haps_file,
            "-s", self.sample_file,
            "-m", self.map_file,
            "-N", str(Ne),
            "-o", str(self.output_prefix),
            "-t", str(self.threads),
        ]
        self._run_command(cmd, "RelateParallel")

    # ==========================================================
    # CONVERT TO tskit TREE SEQUENCE
    # ==========================================================

    def _convert_to_trees(self):
        """
        Convert Relate output to .trees format (tskit-compatible).
        """

        cmd = [
            str(self.convert_bin),
            "--mode",
            "ConvertToTreeSequence",
            "--haps",
            self.haps_file,
            "--sample",
            self.sample_file,
            "--anc",
            f"{self.output_prefix}.anc",
            "--mut",
            f"{self.output_prefix}.mut",
            "--poplabels",
            f"{self.output_prefix}.poplabels",
            "--out",
            self.output_prefix,
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
        total_cM = self.sequence_length * self.recombination_rate * 100.0

        with open(self.map_file, "w") as f:
            f.write(f"{self.chromosome}\tstart\t0.0\t1\n")
            f.write(
                f"{self.chromosome}\tend\t{total_cM:.8f}\t{self.sequence_length}\n"
            )

    # ==========================================================
    # UTILITIES
    # ==========================================================

    def _check_binaries(self):
        if not self.relate_bin.exists():
            raise FileNotFoundError("Relate binary missing")
        if not self.convert_bin.exists():
            raise FileNotFoundError("RelateFileFormats missing")
        if self.mode == "parallel" and not self.parallel_script.exists():
            raise FileNotFoundError("RelateParallel missing")

    def _run_command(self, cmd, label):
        print(f"\n[{label}]")
        print(" ".join(map(str, cmd)))

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(result.stderr)
            raise RuntimeError(f"{label} failed")

        print("done")