import os
import shutil
import subprocess
import numpy as np
import tsinfer
import tsdate
import tskit


class InferTskit:
    def __init__(
        self,
        vcf_file,
        output_prefix,
        output_dir,
        infer_trees=True,
        remove_zarr=True,
        ancestral_mode="from_ts",   
        ts_file=None,
        check_ancestral_state=True,
        num_sites=None,
        apply_tsdate=True,
        tsdate_time_unit="generations",
        tsdate_mutation_rate=1e-8
    ):
        self.vcf_file = vcf_file
        self.ts_file = ts_file
        self.output_prefix = output_prefix
        self.output_dir = output_dir
        self.infer_trees = infer_trees
        self.remove_zarr = remove_zarr
        self.ancestral_mode = ancestral_mode
        #only for from_ts
        self.check_ancestral_state=check_ancestral_state

        self.num_sites = num_sites

        self.ts_inferred = os.path.join(
            output_dir,
            f"{output_prefix}_tsinfer.ts"
        )

        self.zarr_file = os.path.splitext(self.vcf_file)[0] + ".zarr"

        self.apply_tsdate = apply_tsdate
        self.tsdate_time_unit = tsdate_time_unit
        self.tsdate_mutation_rate = tsdate_mutation_rate

    # ==========================================================
    # MAIN
    # ==========================================================

    def run(self):

        if not self.infer_trees:
            return

        print("\n[tsinfer inference]")

        # ------------------------------------------------------
        # Convert VCF to zarr
        # ------------------------------------------------------
        cmd = [
            "vcf2zarr",
            "convert",
            "--force",
            self.vcf_file,
            self.zarr_file,
        ]
        subprocess.run(cmd, check=True)

        print("zarr conversion finished")

        # ------------------------------------------------------
        # Prepare VariantData with ancestral handling
        # ------------------------------------------------------
        vdata = self._prepare_variant_data()

        # ------------------------------------------------------
        # Run tsinfer
        # ------------------------------------------------------
        print("start infer on zarr vdata")
        inferred_ts = tsinfer.infer(vdata)
        print("infer finished")

        if self.apply_tsdate:
            inferred_ts = tsdate.preprocess_ts(
                inferred_ts,
                erase_flanks=False,
                split_disjoint = False
            )
            inferred_ts = tsdate.date(inferred_ts, mutation_rate=self.tsdate_mutation_rate, time_units=self.tsdate_time_unit)

        inferred_ts.dump(self.ts_inferred)
        print("inferred tree dumped:", self.ts_inferred)

        # Cleanup
        if self.remove_zarr:
            if os.path.isfile(self.zarr_file):
                os.remove(self.zarr_file)
            elif os.path.isdir(self.zarr_file):
                shutil.rmtree(self.zarr_file)


    def _prepare_variant_data(self):

        # --------------------------
        # from ts
        if self.ancestral_mode == "from_ts":

            if self.ts_file is None:
                raise ValueError(
                    "ts_file must be provided when ancestral_mode='from_ts'"
                )

            ts_true = tskit.load(self.ts_file)


            ancestral_states = np.zeros(ts_true.num_sites, dtype=object)

            if self.check_ancestral_state:
                for i, site in enumerate(ts_true.sites()):
                    
                    ancestral_states[i] = site.ancestral_state


            return tsinfer.VariantData(self.zarr_file, ancestral_states)

  

        #currently, either num_sites (can be obtained from vcf) or orgiginal ts-file has to be provided
        if not self.num_sites:

            ts_true = tskit.load(self.ts_file)
            ancestral_states = np.zeros(ts_true.num_sites, dtype=object)
        else:
            ancestral_states = np.zeros(self.num_sites, dtype=object)


        vdata = tsinfer.VariantData(self.zarr_file, np.array(ancestral_states))

        if self.ancestral_mode == "none":
            return vdata

        # --------------------------
        # zeros
        if self.ancestral_mode == "zeros":
            ancestral_states = np.zeros(vdata.num_sites, dtype=object)
            return tsinfer.VariantData(self.zarr_file, ancestral_states)

        # --------------------------
        # random - only for binary
        # --------------------------
        elif self.ancestral_mode == "random":
            ancestral_states = np.random.randint(
                0, 2, size=vdata.num_sites
            ).astype(object)
            return tsinfer.VariantData(self.zarr_file, ancestral_states)

        # --------------------------
        # major allele
        # --------------------------
        elif self.ancestral_mode == "major":
            ancestral_states = np.zeros(vdata.num_sites, dtype=object)

            for i, variant in enumerate(vdata.variants()):
                counts = np.bincount(variant.genotypes)
                ancestral_states[i] = np.argmax(counts)

            return tsinfer.VariantData(self.zarr_file, ancestral_states)
        

        raise ValueError(
                f"Unknown ancestral_mode: {self.ancestral_mode}"
            )

