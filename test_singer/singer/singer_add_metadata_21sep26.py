#cd /maps/projects/ilab/people/xcj768/arctic_wolf/args/singer/output/38
#conda activate tskit-env
# inputs: aw.poplabels (metadata), eg trees/aw_n82_38_0.trees
#---------------------------------------------------------------------------------------------------------------------
import tskit
import json
import argparse
#---------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="input params")
parser.add_argument("--ts_file", type=str, help="input ts file after converting singer to tskit format")
parser.add_argument("--poplabels_file", type=str, help="poplabels file")
parser.add_argument("--out_file", type=str, help="out ts file with metadata")
#-------------------------------------------------------------------------------------------------
args = parser.parse_args()
ts_file = args.ts_file
poplabels_file = args.poplabels_file
out_file = args.out_file
#---------------------------------------------------------------------------------------------------------------------
ts = tskit.load(ts_file)
# has 0 inds atm - print(ts.num_individuals)
tables = ts.dump_tables()
# metadata has 4 cols: sample    population group  sex
metadatain = pd.read_csv(poplabels_file, sep=' ')
#---------------------------------------------------------------------------------------------------------------------
tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
# assign number to each pop
pop_name_to_id = {}
for pop_name in metadatain["population"].unique():
    pop_id = tables.populations.add_row(metadata={"name": pop_name})
    pop_name_to_id[pop_name] = pop_id
# update the tables 
tables.individuals.clear()
for row in metadatain.itertuples():
    tables.individuals.add_row(metadata={"name": row.sample})
# nodes to indivduals (1 ind, 2 nodes), assumes diploid 
for i in range(0, ts.num_samples, 2):
    ind_id = i // 2
    pop_name = metadatain.loc[ind_id, "population"]
    target_pop_id = pop_name_to_id[pop_name]
    # node 1 of this ind
    tables.nodes[i] = tables.nodes[i].replace(individual=ind_id, population=target_pop_id)
    # node 2 of this ind
    tables.nodes[i+1] = tables.nodes[i+1].replace(individual=ind_id, population=target_pop_id)
# output this
ts_wmeta = tables.tree_sequence()
ts_wmeta.dump(out_file)
#---------------------------------------------------------------------------------------------------------------------