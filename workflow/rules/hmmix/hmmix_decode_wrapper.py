import os
import sys
import subprocess

HEADER = "chrom\tstart\tend\tlength\tstate\tmean_prob\tsnps\n"


def write_empty(out_file):
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, "w") as f:
        f.write(HEADER)


def main():
    import argparse

    p = argparse.ArgumentParser()

    p.add_argument("--obs")
    p.add_argument("--mutrates")
    p.add_argument("--param")
    p.add_argument("--out")
    p.add_argument("--alpha")
    p.add_argument("--weights", default="")

    #for the empty file we have to manually add the suffix
    suffix_to_add = ".diploid.txt"

    args = p.parse_args()

    required = [args.obs, args.mutrates, args.param, args.alpha]

    # check all inputs exist
    if any(not os.path.exists(x) for x in required):
        print("WARNING: missing inputs, writing empty output", file=sys.stderr)
        empty_file_name = args.out + suffix_to_add
        write_empty(empty_file_name)
        return

    # read alpha
    with open(args.alpha) as f:
        alpha = f.read().strip()

    cmd = [
        "hmmix",
        "decode",
        "-obs", args.obs,
        "-mutrates", args.mutrates,
        "-param", args.param,
        "-out", args.out,
        "-hybrid", alpha
    ]

    if args.weights:
        cmd += ["-weights", args.weights]

    subprocess.run(cmd, check=False)


if __name__ == "__main__":
    main()