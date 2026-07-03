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


def build_trace_indices(ref_file, tgt_file, ploidy=2):
    ref = load_ids(ref_file)
    tgt = load_ids(tgt_file)

    ref_idx = expand_individuals(ref, ploidy, start_index=0)
    tgt_idx = expand_individuals(tgt, ploidy, start_index=len(ref_idx))

    return ref_idx, tgt_idx


if __name__ == "__main__":
    ref_file = sys.argv[1]
    tgt_file = sys.argv[2]
    ploidy = int(sys.argv[3]) if len(sys.argv) > 3 else 2

    _, tgt_idx = build_trace_indices(ref_file, tgt_file, ploidy)

    print(",".join(map(str, tgt_idx)))