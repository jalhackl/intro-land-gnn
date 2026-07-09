#-----------------------------------------------------------------------------------------------------------------------
# code to add metadata to relate trees from -
# https://github.com/leospeidel/relate_lib/issues/7
#-----------------------------------------------------------------------------------------------------------------------
import tskit
import tszip
import json
import argparse
#-----------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="input params")
parser.add_argument("--ts_file", type=str, help="input ts file from relate_lib conversion")
parser.add_argument("--poplabels_file", type=str, help="poplabels file")
parser.add_argument("--out_prefix", type=str, help="out prefix")
#-------------------------------------------------------------------------------------------------
args = parser.parse_args()
#---------------------------------------------------------------------------------------------------------------------
def dict2byte_string(d):
    return json.dumps(d).encode("utf-8")
def byte_string2dict(b):
    return json.loads(b.decode("utf-8"))
#-----------------------------------------------------------------------------------------------------------------------
def add_pop_indiv_info(args):
    ts_file = args.ts_file
    poplabels_file = args.poplabels_file
    out_prefix = args.out_prefix
    print(f'ts_file: {ts_file}')
    print(f'poplabels_file: {poplabels_file}')
    print(f'out_prefix: {out_prefix}')
    print("Loading tree sequence from %s"%(ts_file))
    ts = tszip.load(ts_file)
    tables = ts.dump_tables()
    ## load .poplabels
    print(f'Loading population and individual info from {poplabels_file}')
    uniq_pops = []
    #
    metadata_column = []
    with open(poplabels_file, 'r') as ifl:
        header = ifl.readline().strip().split()
        assert header == ['sample', 'population', 'group', 'sex']
        #
        for l in ifl:
            line = l.strip().split()
            metadata = {}
            for k, v in zip(header, line):
                metadata[k] = v
                #
            metadata_column.append(dict2byte_string(metadata))
            uniq_pops.append(line[1])
            #
    uniq_pops = sorted(set(uniq_pops))
    print(f'- Number of individuals: {len(metadata_column)}')
    print(f'- Number of unique populations: {len(uniq_pops)}')
    #
    ## set individuals table
    print("Setting individuals table")
    # basic_schema = tskit.MetadataSchema({'codec': 'json'})
    basic_schema = tskit.MetadataSchema(None)
    tables.individuals.metadata_schema = basic_schema
    #
    flags = [0 for x in range(len(metadata_column))]
    #
    encoded_metadata_column = [
        tables.individuals.metadata_schema.validate_and_encode_row(r) for r in metadata_column
    ]
    md, md_offset = tskit.pack_bytes(encoded_metadata_column)
    tables.individuals.set_columns(flags=flags, metadata=md, metadata_offset=md_offset)
    #
    ## set populations table
    print("Setting populations table")
    population_metadata_column = []
    #
    for pop in uniq_pops:
        population_metadata_column.append(dict2byte_string({'name': pop}))
        #
    tables.populations.metadata_schema = basic_schema
    #
    encoded_pop_metadata_column = [
        tables.populations.metadata_schema.validate_and_encode_row(r) for r in population_metadata_column
    ]
    md, md_offset = tskit.pack_bytes(encoded_pop_metadata_column)
    tables.populations.set_columns(metadata=md, metadata_offset=md_offset)
    #
    ## set individual id of sample nodes in node table
    print("Setting individual id of sample nodes in node table")
    num_individuals = len(metadata_column)
    new_individual_of_sample_nodes = []
    for i in range(num_individuals):
        new_individual_of_sample_nodes.extend([i, i])
        #
    new_individual = tables.nodes.individual
    new_individual[:len(new_individual_of_sample_nodes)] = new_individual_of_sample_nodes
    #
    tables.nodes.individual = new_individual
    #
    ## set population of sample nodes in node table
    print("Setting population id of sample nodes in node table")
    num_indivs = len(tables.individuals)
    new_population_for_sample_nodes = []
    for idx in range(num_indivs):
        indiv = tables.individuals[idx]
        pop = byte_string2dict(indiv.metadata)['population']
        pop_idx = uniq_pops.index(pop)
        new_population_for_sample_nodes.extend([pop_idx, pop_idx])
        #
    new_population = tables.nodes.population
    new_population[:len(new_population_for_sample_nodes)] = new_population_for_sample_nodes
    #
    tables.nodes.population = new_population
    #
    ## convert to new tree sequence
    new_ts = tables.tree_sequence()
    #
    out_trees_file = f'{out_prefix}.trees.tsz'
    tszip.compress(new_ts, out_trees_file)
    #
    print(f'Done. out prefix: {out_prefix}')
    print(f'Output tree sequence: {out_trees_file}')
#-----------------------------------------------------------------------------------------------------------------------
add_pop_indiv_info(args)
#-----------------------------------------------------------------------------------------------------------------------