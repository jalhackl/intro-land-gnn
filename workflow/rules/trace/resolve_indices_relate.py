import sys

def load_ids(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def expand_individuals(individuals, ploidy, start_index=0):
    """
    Returns haplotype indices for individuals.

    Each individual i maps to:
        start_index + i*ploidy ... start_index + i*ploidy + (ploidy-1)
    """
    idx = []
    current = start_index

    for _ in individuals:
        idx.extend(range(current, current + ploidy))
        current += ploidy

    return idx


def build_trace_indices_simple(ref_file, tgt_file, ploidy=2):
    ref = load_ids(ref_file)
    tgt = load_ids(tgt_file)

    ref_idx = expand_individuals(ref, ploidy, start_index=0)
    tgt_idx = expand_individuals(tgt, ploidy, start_index=len(ref_idx))

    return ref_idx, tgt_idx




def build_trace_indices(ref_file, tgt_file, ploidy=2, haplotype_start_index=0):
    ref = load_ids(ref_file)
    tgt = load_ids(tgt_file)

    ref_idx = expand_individuals(ref, ploidy, start_index=0)

    tgt_id_map = dict()
    id_tgt_map = dict()

    if ploidy > 1:
        ref_pl = []
        tgt_pl = []
        for ref_ind in ref:
            current = haplotype_start_index
            for h in range(ploidy):
                name = f"{ref_ind}_{h+1}"
                current += 1
                ref_pl.append(name)
        ref = ref_pl


    tgt_idx = expand_individuals(tgt, ploidy, start_index=len(ref_idx))



    curr_tgt_idx = tgt_idx[0]
    #new names
    if ploidy > 1:

        tgt_pl = []

        for tgt_ind in tgt:
            current = haplotype_start_index
            for h in range(ploidy):
                name = f"{tgt_ind}_{h+1}"
                current += 1
                tgt_pl.append(tgt)
                tgt_id_map[name] = str(curr_tgt_idx)
                id_tgt_map[str(curr_tgt_idx)] = name

                curr_tgt_idx += 1

    else:
        for tgt_ind in tgt:
            name = tgt_ind
            tgt_id_map[name] = str(curr_tgt_idx)
            id_tgt_map[str(curr_tgt_idx)] = name

            curr_tgt_idx += 1



    return ref_idx, tgt_idx, id_tgt_map, tgt_id_map

if __name__ == "__main__":
    ref_file = sys.argv[1]
    tgt_file = sys.argv[2]
    ploidy = int(sys.argv[3]) if len(sys.argv) > 3 else 2

    _, tgt_idx, tgt_id_map, tgt_name_to_index = build_trace_indices(ref_file, tgt_file, ploidy)

    print(",".join(map(str, tgt_idx)))