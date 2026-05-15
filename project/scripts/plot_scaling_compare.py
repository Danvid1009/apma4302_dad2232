"""Overlay speedup (and optional efficiency) curves from multiple scaling CSVs."""

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Compare speedup curves from several scaling CSVs.")
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Each entry: path.csv:Label e.g. output/a.csv:European",
    )
    parser.add_argument("--out", type=str, default="output/scaling_compare_speedup.png")
    parser.add_argument("--title", type=str, default="Scaling speedup comparison")
    return parser.parse_args()


def load_series(csv_path: pathlib.Path):
    ranks, speedup, efficiency = [], [], []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ranks.append(int(row["ranks"]))
            speedup.append(float(row["speedup"]))
            efficiency.append(float(row["efficiency"]))
    return ranks, speedup, efficiency


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure()
    max_rank = 1
    for spec in args.inputs:
        if ":" in spec:
            path_s, label = spec.split(":", 1)
        else:
            path_s, label = spec, spec
        ranks, speedup, _ = load_series(pathlib.Path(path_s))
        max_rank = max(max_rank, max(ranks))
        plt.plot(ranks, speedup, marker="o", label=label)

    plt.plot([1, max_rank], [1, max_rank], linestyle="--", color="gray", label="Ideal")
    plt.xlabel("MPI ranks")
    plt.ylabel("Speedup")
    plt.title(args.title)
    plt.legend()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
