# issues converting output of singer to tskit format

#-----------------------------------------------------------------------------------------------------------------------

(tskit-env) [xcj768@mjolnirgate01fl arctic_wolves]$ grep 'error' /projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/parallelsinger_43194318-38.log
(tskit-env) [xcj768@mjolnirgate01fl arctic_wolves]$ grep 'Error' /projects/psg/people/xcj768/arctic_wolves/impute/logs/final_dataset/parallelsinger_43194318-38.log

# parallel singer looks ok
#-----------------------------------------------------------------------------------------------------------------------



conda activate tskit-env

module load samtools/1.21 tabix/1.11 
module load  snakemake/9.9.0 
module load gsl/2.5 perl bcftools/1.21  



scp -r  \
~/Library/CloudStorage/OneDrive-UniversityofCopenhagen/CPH_scripts/cph_arctic_wolf/converttotskit.arr \
xcj768@mjolnirgate.unicph.domain:/projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg

cd /projects/psg/people/xcj768/arctic_wolves
sbatch /projects/psg/people/xcj768/arctic_wolves/scripts/downstream/final_dataset/arg/converttotskit.arr
# Submitted batch job 43220267 Tue Jan 20 14:02:00 CET 2026

ID=38

module load parallel/20241222 

ID=38
input="/projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/${ID}/aw_n82_${ID}"
 outp="/projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/${ID}/trees/aw_n82_${ID}"
/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit \
>  -input ${input} \
>  -output ${outp} \
>  -start 0 \
>  -end 99
Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 98, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 95, in main
    write_trees(args.input, args.output, args.start, args.end, args.step)
    ~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 64, in write_trees
    ts = read_ARG(node_file, branch_file, mutation_file)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 52, in read_ARG
    tables = read_ts(node_file, branch_file)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 10, in read_ts
    node_time = np.loadtxt(node_file)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1384, in loadtxt
    arr = _read(fname, dtype=dtype, comment=comment, delimiter=delimiter,
                converters=converters, skiplines=skiprows, usecols=usecols,
                unpack=unpack, ndmin=ndmin, encoding=encoding,
                max_rows=max_rows, quote=quotechar)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1011, in _read
    fh = np.lib._datasource.open(fname, 'rt', encoding=encoding)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 192, in open
    return ds.open(path, mode, encoding=encoding, newline=newline)
           ~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 529, in open
    raise FileNotFoundError(f"{path} not found.")
FileNotFoundError: /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_nodes_0.txt not found.

(tskit-env) [xcj768@mjolnirgate01fl arctic_wolves]$ ls -lhtr /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38_nodes_0.txt
ls: cannot access '/projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38_nodes_0.txt': No such file or directory


(tskit-env) [xcj768@mjolnirgate01fl arctic_wolves]$ ls -lhtr  /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/*nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 199K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_3_4_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 266K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_0_1_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 264K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_5_6_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 265K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_1_2_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 281K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_2_3_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 273K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_14_15_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 250K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_16_17_start_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 191K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_3_4_nodes_0.txt
-rw-r--r--. 1 xcj768 comp-prj-psg 245K Jan 19 22:34 /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_9_10_start_nodes_0.txt

# whether this needs to be diff

prefix_nodes_0.txt but prefix_0_1_nodes_0.txt

#-----------------------------------------------------------------------------------------------------------------------
https://github.com/popgenmethods/SINGER/issues/21

@zhangming-m there is a "convert_long_ARG.py" in the repo which you might want to use. So the usage is like 
"python convert_long_ARG.py vcf_prefix output_prefix iteration_index". 

If this does not fix your problem, you can send me the files and I can convert it for you. 
Just find all the vcf_prefix*_nodes_0.txt, vcf_prefix*_branches_0.txt, vcf_prefix*_muts_0.txt, vcf_prefix*_recombs_0.txt and send them to me. I can work from there.


less /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py
#-----------------------------------------------------------------------------------------------------------------------

# try to run convert_long_ARG.py
"python convert_long_ARG.py vcf_prefix output_prefix iteration_index". 


(tskit-env) [xcj768@mjolnirgate01fl arctic_wolves]$ head /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38.index
0	8092
1000000	8091904
2000000	16479604


python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
args/singer/output/test/38/aw_n82_38 \
args/singer/output/test/38/trees \
0

usage: convert_long_ARG.py [-h] -vcf VCF -output OUTPUT -iteration ITERATION
convert_long_ARG.py: error: the following arguments are required: -vcf, -output, -iteration

cd args/singer/output/test/38/
python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
-vcf aw_n82_38 \
-output trees/aw_n82_38_0 \
-iteration 0

Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 118, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 111, in main
    node_files, branch_files, mutation_files, block_coordinates = generate_file_lists(args.vcf, args.output, args.iteration)
                                                                  ~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 80, in generate_file_lists
    with open(f"{vcf_prefix}.index", 'r') as f:
         ~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FileNotFoundError: [Errno 2] No such file or directory: 'aw_n82_38.index'

# add soft link 
ln -s /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38.index .

python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
-vcf aw_n82_38 \
-output trees/aw_n82_38_0 \
-iteration 0
Processing segment 0
Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 118, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 114, in main
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 23, in read_long_ARG
    node_time = np.loadtxt(node_file)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1384, in loadtxt
    arr = _read(fname, dtype=dtype, comment=comment, delimiter=delimiter,
                converters=converters, skiplines=skiprows, usecols=usecols,
                unpack=unpack, ndmin=ndmin, encoding=encoding,
                max_rows=max_rows, quote=quotechar)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1011, in _read
    fh = np.lib._datasource.open(fname, 'rt', encoding=encoding)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 192, in open
    return ds.open(path, mode, encoding=encoding, newline=newline)
           ~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 529, in open
    raise FileNotFoundError(f"{path} not found.")
FileNotFoundError: trees/aw_n82_38_0_0_1_nodes_0.txt not found.
#-----------------------------------------------------------------------------------------------------------------------
python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py -h

usage: convert_long_ARG.py [-h] -vcf VCF -output OUTPUT -iteration ITERATION

Generate tskit format for a long ARG.

options:
  -h, --help            show this help message and exit
  -vcf VCF              VCF file prefix
  -output OUTPUT        Output files prefix
  -iteration ITERATION  MCMC iteration for generating filenames]

# confused by this..
#-----------------------------------------------------------------------------------------------------------------------
/projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38.vcf

python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
-vcf /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38 \
-output aw_n82_38 \
-iteration 0



(tskit-env) [xcj768@mjolnirgate01fl 38]$ python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
> -vcf /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38 \
> -output aw_n82_38 \
> -iteration 0
Processing segment 0
Processing segment 1
Processing segment 2
Processing segment 3
Processing segment 4
Processing segment 5
Processing segment 6
Processing segment 7
Processing segment 8
Processing segment 9
Processing segment 10
Processing segment 11
Processing segment 12
Processing segment 13
Processing segment 14
Processing segment 15
Processing segment 16
Processing segment 17
Processing segment 18
Processing segment 19
Processing segment 20
Processing segment 21
Processing segment 22
Processing segment 23
Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 118, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 114, in main
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 74, in read_long_ARG
    ts = tables.tree_sequence()
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/tskit/tables.py", line 3528, in tree_sequence
    return tskit.TreeSequence.load_tables(self)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/tskit/trees.py", line 4295, in load_tables
    ts.load_tables(tables._ll_tables, build_indexes=build_indexes)
    ~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#_tskit.LibraryError: A mutation's parent is not consistent with the topology of the tree. Use compute_mutation_parents to set the parents correctly.(TSK_ERR_BAD_MUTATION_PARENT)


# make github issue 
#-----------------------------------------------------------------------------------------------------------------------


# may need to write github issue * - have writtend
- issue converting to tskit format 

- ran parallel singer per chrom


output as follows
(tskit-env) [xcj768@mjolnirgate01fl 38]$ ls -lhtr  aw_n82_38_0_1*_0.txt
 aw_n82_38_0_1_start_nodes_0.txt
 aw_n82_38_0_1_start_branches_0.txt
 aw_n82_38_0_1_start_recombs_0.txt
 aw_n82_38_0_1_start_muts_0.txt
 aw_n82_38_0_1_nodes_0.txt
 aw_n82_38_0_1_branches_0.txt
 aw_n82_38_0_1_recombs_0.txt
 aw_n82_38_0_1_muts_0.txt

- attempted to convert to tskit format
2 ways

q1) which is the recommended method to convert to tskit format after running parallel singer
using convert_to_tskit or convert_long_ARG.py?

# i guess b/c doesn't expect per window output 
ID=38
input="/projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/${ID}/aw_n82_${ID}"
 outp="/projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/${ID}/trees/aw_n82_${ID}"
/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit \
>  -input ${input} \
>  -output ${outp} \
>  -start 0 \
>  -end 99
Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 98, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 95, in main
    write_trees(args.input, args.output, args.start, args.end, args.step)
    ~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 64, in write_trees
    ts = read_ARG(node_file, branch_file, mutation_file)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 52, in read_ARG
    tables = read_ts(node_file, branch_file)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_to_tskit", line 10, in read_ts
    node_time = np.loadtxt(node_file)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1384, in loadtxt
    arr = _read(fname, dtype=dtype, comment=comment, delimiter=delimiter,
                converters=converters, skiplines=skiprows, usecols=usecols,
                unpack=unpack, ndmin=ndmin, encoding=encoding,
                max_rows=max_rows, quote=quotechar)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_npyio_impl.py", line 1011, in _read
    fh = np.lib._datasource.open(fname, 'rt', encoding=encoding)
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 192, in open
    return ds.open(path, mode, encoding=encoding, newline=newline)
           ~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/numpy/lib/_datasource.py", line 529, in open
    raise FileNotFoundError(f"{path} not found.")
FileNotFoundError: /projects/psg/people/xcj768/arctic_wolves/args/singer/output/test/38/aw_n82_38_nodes_0.txt not found.

=> following https://github.com/popgenmethods/SINGER/issues/21 I tried convert_long_ARG.py

q2) how to fix the error seen in convert script


python /projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py \
> -vcf /projects/psg/people/xcj768/arctic_wolves/args/singer/input/aw_n82_38 \
> -output aw_n82_38 \
> -iteration 0
Processing segment 0
Processing segment 1
Processing segment 2
Processing segment 3
Processing segment 4
Processing segment 5
Processing segment 6
Processing segment 7
Processing segment 8
Processing segment 9
Processing segment 10
Processing segment 11
Processing segment 12
Processing segment 13
Processing segment 14
Processing segment 15
Processing segment 16
Processing segment 17
Processing segment 18
Processing segment 19
Processing segment 20
Processing segment 21
Processing segment 22
Processing segment 23
Traceback (most recent call last):
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 118, in <module>
    main()
    ~~~~^^
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 114, in main
    ts = read_long_ARG(node_files, branch_files, mutation_files, block_coordinates)
  File "/projects/psg/people/xcj768/resources/singer_v0.1.8/singer-0.1.8-beta-linux-x86_64/convert_long_ARG.py", line 74, in read_long_ARG
    ts = tables.tree_sequence()
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/tskit/tables.py", line 3528, in tree_sequence
    return tskit.TreeSequence.load_tables(self)
           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^
  File "/home/xcj768/miniforge3/envs/tskit-env/lib/python3.14/site-packages/tskit/trees.py", line 4295, in load_tables
    ts.load_tables(tables._ll_tables, build_indexes=build_indexes)
    ~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
_tskit.LibraryError: A mutation's parent is not consistent with the topology of the tree. Use compute_mutation_parents to set the parents correctly.(TSK_ERR_BAD_MUTATION_PARENT)




#-----------------------------------------------------------------------------------------------------------------------
