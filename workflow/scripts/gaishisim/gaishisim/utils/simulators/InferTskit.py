import os
import shutil
import subprocess
import numpy as np
import tsinfer
import tsdate
import tskit
import json
import zarr
import pandas as pd

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
        fasta_file=None,
        check_ancestral_state=True,
        num_sites=None,
        apply_tsdate=True,
        tsdate_time_unit="generations",
        tsdate_mutation_rate=1.2e-8,
        
        #zarr_variants_chunk_size=None,
        #zarr_samples_chunk_size=None,
        #worker_processes=None,
        zarr_variants_chunk_size=100000,
        zarr_samples_chunk_size=1000,
        worker_processes=4,
        use_pyfaidx=False,
        #poplabels_file="subset_individuals_relate_formatting.txt",
        poplabels_file=None,
        erase_flanks=True,
        split_disjoint=True,
        check_existing_zarr=True
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
        
        self.ts_inferred_before_tsdate = os.path.join(
            output_dir,
            f"{output_prefix}_tsinfer_before_tsdate.ts"
        )


        self.zarr_file = os.path.splitext(self.vcf_file)[0] + ".zarr"
        self.poplabels_file = poplabels_file

        self.apply_tsdate = apply_tsdate
        self.tsdate_time_unit = tsdate_time_unit
        self.tsdate_mutation_rate = tsdate_mutation_rate
        
        self.fasta_file = fasta_file
        
        self.zarr_variants_chunk_size = zarr_variants_chunk_size
        self.zarr_samples_chunk_size = zarr_samples_chunk_size
        self.worker_processes = worker_processes
        self.use_pyfaidx = use_pyfaidx
        
        self.erase_flanks = erase_flanks
        self.split_disjoint = split_disjoint
        self.check_existing_zarr = check_existing_zarr

    # ==========================================================
    # MAIN
    # ==========================================================

    def run(self):

        if not self.infer_trees:
            return

        print("\n[tsinfer inference]")


        if self.check_existing_zarr and os.path.exists(self.zarr_file):
            print(f"Zarr file {self.zarr_file} already exists. Skipping conversion.") 
            
        else:  
            # ------------------------------------------------------
            # Convert VCF to zarr
            # ------------------------------------------------------
            '''
            cmd = [
                "vcf2zarr",
                "convert",
                "--force",
                self.vcf_file,
                self.zarr_file,
            ]
            '''
            
            cmd = [
                "vcf2zarr",
                "convert",
                "--force",
            ]

            if self.zarr_variants_chunk_size is not None:
                cmd.extend([
                    "--variants-chunk-size",
                    str(self.zarr_variants_chunk_size)
                ])

            if self.zarr_samples_chunk_size is not None:
                cmd.extend([
                    "--samples-chunk-size",
                    str(self.zarr_samples_chunk_size)
                ])

            if self.worker_processes is not None:
                cmd.extend([
                    "--worker-processes",
                    str(self.worker_processes)
                ])

            cmd.extend([
                self.vcf_file,
                self.zarr_file,
            ])
                    
            
            subprocess.run(cmd, check=True)

            print("zarr conversion finished")
            
            if self.poplabels_file is not None:
                self._add_population_labels()

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


        inferred_ts.dump(self.ts_inferred_before_tsdate)


        if self.apply_tsdate:
            inferred_ts = tsdate.preprocess_ts(
                inferred_ts,
                erase_flanks=self.erase_flanks,
                split_disjoint=self.split_disjoint
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
        
        # from FASTA
        if self.ancestral_mode == "fasta":

            if self.fasta_file is None:
                raise ValueError(
                    "fasta_file must be provided when ancestral_mode='fasta'"
                )

            if self.use_pyfaidx:
                import pyfaidx
                import zarr

                vcf_zarr = zarr.open(self.zarr_file)

                reader = pyfaidx.Fasta(self.fasta_file)

                # assuming ancestral chromosome FASTA has only one sequence
                seqname = reader.keys()[0]

                ancestral_sequence = str(reader[seqname]).upper()
                
                # Replace anything not A/C/G/T with N
                ancestral_sequence = "".join(
                    base if base in {"A", "C", "G", "T"} else "N"
                    for base in ancestral_sequence
                )

                # Variant positions are one-based, so prepend dummy character
                ancestral_sequence = "X" + ancestral_sequence
                
                # Convert to numpy array of single-character strings
                ancestral_sequence = np.array(list(ancestral_sequence), dtype="U1")

                tsinfer.add_ancestral_state_array(
                    vcf_zarr,
                    ancestral_sequence
                )

                return tsinfer.VariantData(
                    self.zarr_file,
                    ancestral_state="ancestral_state"
                )

            else:
                import pysam
                
                vcf = pysam.VariantFile(self.vcf_file)
                fasta = pysam.FastaFile(self.fasta_file)

                # assuming ancestral chromosome FASTA has only one sequence
                seqname = fasta.references[0]

                #ancestral_sequence = fasta.fetch(seqname)
                
                
                #ancestral_states = [
                #    fasta.fetch(seqname, record.pos - 1, record.pos)
                #    for record in vcf
                #]
                
                
                ancestral_states = []

                for record in vcf:
                    base = fasta.fetch(seqname, record.pos - 1, record.pos).upper()
                    ancestral_states.append(base if base in {"A", "C", "G", "T"} else "N")

                ancestral_states = np.array(ancestral_states, dtype="U1")
                
                
                '''
                ancestral_states = np.zeros(
                    self.num_sites,
                    dtype=object
                )

                vdata = tsinfer.VariantData(self.zarr_file, ancestral_states)

                for i, variant in enumerate(vdata.variants()):
                    ancestral_states[i] = ancestral_sequence[variant.position - 1]
                '''
                return tsinfer.VariantData(
                    self.zarr_file,
                    ancestral_states
                )
        
            
        # from VCF AA field
        if self.ancestral_mode == "from_vcf_aa":

            vdata = tsinfer.VariantData(
                self.zarr_file,
                np.zeros(self.num_sites, dtype=object)
            )

            ancestral_states = np.zeros(
                vdata.num_sites,
                dtype=object
            )

            vcf = pysam.VariantFile(self.vcf_file)

            for i, record in enumerate(vcf):

                if "AA" not in record.info:
                    raise ValueError(
                        f"VCF record {record.chrom}:{record.pos} has no AA INFO field"
                    )

                aa = record.info["AA"]

                if isinstance(aa, tuple):
                    aa = aa[0]

                ancestral_states[i] = aa

            return tsinfer.VariantData(
                self.zarr_file,
                ancestral_states
    )

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


    def _add_population_labels(self, population_identifier="group"):
        """
        Add population assignments and metadata to the tsinfer zarr file.
        """
        
        '''
        if treat_group_as_population:
            population_identifier = "group"
        else:
            population_identifier = "population"
        '''

        if self.poplabels_file is None:
            return

        print("Adding population labels to zarr")

        # Read whitespace-separated file
        for encoding in ("utf-8", "iso-8859-1", "cp1252"):
            try:
                labels = pd.read_csv(
                    self.poplabels_file,
                    sep=r"\s+",
                    encoding=encoding,
                    dtype=str,
                )
                break
            except UnicodeDecodeError:
                pass
        else:
            raise ValueError(
                f"Could not decode {self.poplabels_file}"
            )

        required = {"sample", population_identifier}

        if not required.issubset(labels.columns):
            raise ValueError(
                f"poplabels file must contain columns: {required}"
            )

        # Open zarr
        vcf_zarr = zarr.open(self.zarr_file, mode="a")

        # --------------------------------------------------
        # Add population metadata
        # --------------------------------------------------

        populations = labels[population_identifier].unique().tolist()

        schema = json.dumps(
            tskit.MetadataSchema.permissive_json().schema
        ).encode()

        zarr.save(
            f"{self.zarr_file}/populations_metadata_schema",
            schema
        )

        population_metadata = [
            json.dumps(
                {
                    "name": pop
                }
            ).encode()
            for pop in populations
        ]

        zarr.save(
            f"{self.zarr_file}/populations_metadata",
            population_metadata
        )


        # --------------------------------------------------
        # Match VCF samples to population
        # --------------------------------------------------

        sample_ids = [
            s.decode() if isinstance(s, bytes) else s
            for s in vcf_zarr["sample_id"]
        ]

        individuals_population = np.full(
            len(sample_ids),
            tskit.NULL,
            dtype=np.int32
        )


        sample_to_population = dict(
            zip(
                labels["sample"],
                labels[population_identifier]
            )
        )


        for i, sample in enumerate(sample_ids):

            if sample not in sample_to_population:
                raise ValueError(
                    f"Sample {sample} missing in poplabels file"
                )

            pop = sample_to_population[sample]

            individuals_population[i] = populations.index(pop)


        zarr.save(
            f"{self.zarr_file}/individuals_population",
            individuals_population
        )


        # --------------------------------------------------
        # Optional individual metadata
        # --------------------------------------------------

        metadata = []

        label_dict = labels.set_index("sample").to_dict(
            orient="index"
        )

        for sample in sample_ids:

            if sample in label_dict:
                metadata.append(
                    json.dumps(
                        label_dict[sample]
                    ).encode()
                )
            else:
                metadata.append(
                    json.dumps({}).encode()
                )

        zarr.save(
            f"{self.zarr_file}/individuals_metadata",
            metadata
        )

        print(
            f"Added {len(populations)} populations"
        )
