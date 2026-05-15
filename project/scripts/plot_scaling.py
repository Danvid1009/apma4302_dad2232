import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Plot speedup/efficiency from scaling CSV.")
    parser.add_argument("--infile", type=str, default="output/scaling_results.csv")
    parser.add_argument("--out-prefix", type=str, default="output/scaling")
    return parser.parse_args()


def main():
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    out_prefix = pathlib.Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ranks, speedup, efficiency = [], [], []
    with in_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ranks.append(int(row["ranks"]))
            speedup.append(float(row["speedup"]))
            efficiency.append(float(row["efficiency"]))

    plt.figure()
    plt.plot(ranks, speedup, marker="o", label="Measured")
    plt.plot(ranks, ranks, linestyle="--", label="Ideal")
    plt.xlabel("MPI ranks")
    plt.ylabel("Speedup")
    plt.title("Strong/Weak Scaling Speedup")
    plt.legend()
    speedup_path = f"{out_prefix}_speedup.png"
    plt.savefig(speedup_path, dpi=150, bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(ranks, efficiency, marker="o")
    plt.xlabel("MPI ranks")
    plt.ylabel("Parallel efficiency")
    plt.title("Scaling Efficiency")
    plt.ylim(0, 1.05)
    efficiency_path = f"{out_prefix}_efficiency.png"
    plt.savefig(efficiency_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"wrote {speedup_path}")
    print(f"wrote {efficiency_path}")


if __name__ == "__main__":
    main()
