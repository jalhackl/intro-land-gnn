import numpy as np
import msprime
import yaml



def generate_mutation_rate_map(
    genome_length,
    mean_block_length=100000,
    mode="exponential",
    # Gamma-related params (existing modes)
    mean_rate=1.25e-8,
    # Coefficient of variation as in Cousins
    cv=0.15,
    # shape parameter
    gamma_shape=2.5,
    # only for normal
    normal_sd=1e-9,
    seed=None,
):
    """
    mutation-rate map generator supporting three models:

    1. mode="exponential":
       - block lengths ~ Exponential(mean = mean_block_length)
       - rates         ~ Gamma(mean = mean_rate, cv = cv)

    2. mode="geometric":
       - block lengths ~ Geometric(mean = mean_block_length)
       - rates         ~ Gamma(mean = mean_rate, shape = gamma_shape)

    3. mode="fixed_normal":
       - block lengths = mean_block_length (constant)
       - rates         ~ Normal(mean_rate, normal_sd)
         (rates are truncated at > 0 to avoid negative mutation rates)
    """

    rng = np.random.default_rng(seed)

    positions = [0.0]
    rates = []
    pos = 0.0

    if mode == "exponential":
        shape = 1.0 / (cv**2)
        scale = mean_rate / shape

    elif mode == "geometric":
        shape = gamma_shape
        scale = mean_rate / gamma_shape
        p = 1.0 / mean_block_length

    else:
        raise ValueError("mode must be 'exponential', 'geometric', or 'fixed_normal'")

    # segments

    while pos < genome_length:

        # length of segment
        if mode == "exponential":
            seg_len = rng.exponential(scale=mean_block_length)

        elif mode == "geometric":
            seg_len = rng.geometric(p)

        else:  # normal
            seg_len = mean_block_length

        end = min(pos + seg_len, genome_length)

        # mutation rate
        if mode in ("exponential", "geometric"):
            rate = rng.gamma(shape=shape, scale=scale)
        else:  # normal
            rate = rng.normal(mean_rate, normal_sd)
            while rate < 0:
                rate = rng.normal(mean_rate, normal_sd)

        positions.append(float(end))
        rates.append(float(rate))

        pos = end

    positions = np.asarray(positions)
    rates = np.asarray(rates)

    rate_map = msprime.RateMap(position=positions, rate=rates)

    return rate_map
    # return np.asarray(positions), np.asarray(rates)




def effective_population_size(
    yaml_path,
    mode="all",
    demes_subset=None,
    summary="extreme"   # default changed to "extreme"
):
    """
    Compute effective population size (Ne) from a demes YAML file.

    Parameters
    ----------
    yaml_path : str
        Path to YAML file.

    mode : str
        - "all"     : time-averaged Ne per deme
        - "archaic" : oldest deme(s)
        - "recent"  : most recent deme(s)
        - "subset"  : only selected demes (time-averaged)

    demes_subset : list[str] or None
        Optional filter applied before computation.

    summary : str
        - "mean"    : time-averaged Ne across epochs
        - "extreme" : boundary value
                      (archaic → earliest size,
                       recent → latest size; DEFAULT)
    """

    with open(yaml_path, "r") as f:
        model = yaml.safe_load(f)

    demes = model["demes"]

    # -----------------------------
    # GLOBAL FILTER
    # -----------------------------
    if demes_subset is not None:
        demes = [d for d in demes if d["name"] in demes_subset]

    # -----------------------------
    # Helper: time-averaged Ne
    # -----------------------------
    def time_averaged_ne(deme):
        epochs = deme.get("epochs", [])

        total_time = 0.0
        harmonic_sum = 0.0

        for ep in epochs:
            start = ep.get("start_time", None)
            end = ep.get("end_time", 0.0)

            if start is None:
                continue

            duration = start - end
            if duration <= 0:
                continue

            start_size = ep.get("start_size")
            end_size = ep.get("end_size")

            if end_size is None or start_size == end_size:
                Ne = start_size
            else:
                Ne = 0.5 * (start_size + end_size)

            total_time += duration
            harmonic_sum += duration / Ne

        if harmonic_sum == 0:
            return None

        return total_time / harmonic_sum

    # -----------------------------
    # ALL / SUBSET MODE
    # -----------------------------
    if mode in ("all", "subset"):
        return {d["name"]: time_averaged_ne(d) for d in demes}

    # -----------------------------
    # ARCHAIC MODE
    # -----------------------------
    if mode == "archaic":
        max_start = max(
            ep.get("start_time", 0)
            for d in demes
            for ep in d.get("epochs", [])
            if ep.get("start_time") is not None
        )

        archaic_demes = [
            d for d in demes
            for ep in d.get("epochs", [])
            if ep.get("start_time") == max_start
        ]

        result = {}

        for d in archaic_demes:
            if summary == "mean":
                result[d["name"]] = time_averaged_ne(d)

            elif summary == "extreme":
                # earliest population size (oldest epoch)
                ep = max(d["epochs"], key=lambda e: e.get("start_time", 0))
                result[d["name"]] = ep.get("start_size")

            else:
                raise ValueError("summary must be 'mean' or 'extreme'")

        return result

    # -----------------------------
    # RECENT MODE (DEFAULT EXTREME)
    # -----------------------------
    if mode == "recent":
        min_end = min(
            ep.get("end_time", 0)
            for d in demes
            for ep in d.get("epochs", [])
        )

        recent_demes = [
            d for d in demes
            for ep in d.get("epochs", [])
            if ep.get("end_time", 0) == min_end
        ]

        result = {}

        for d in recent_demes:
            if summary == "mean":
                result[d["name"]] = time_averaged_ne(d)

            elif summary == "extreme":
                # latest population size (boundary at present)
                ep = min(d["epochs"], key=lambda e: e.get("end_time", float("inf")))

                result[d["name"]] = (
                    ep.get("end_size")
                    if ep.get("end_size") is not None
                    else ep.get("start_size")
                )

            else:
                raise ValueError("summary must be 'mean' or 'extreme'")

        return result

    raise ValueError("mode must be one of: all, subset, archaic, recent")