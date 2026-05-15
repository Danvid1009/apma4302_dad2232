"""Plot price and stderr vs barrier from run_barrier_sweep.py CSV."""

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description="Plot barrier sweep CSV.")
    p.add_argument("--infile", type=str, default="output/barrier_sweep.csv")
    p.add_argument("--out-prefix", type=str, default="output/barrier_sweep")
    return p.parse_args()


def main():
    args = parse_args()
    in_path = pathlib.Path(args.infile)
    out_prefix = pathlib.Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    b, price, stderr = [], [], []
    with in_path.open("r", newline="") as f:
        for row in csv.DictReader(f):
            b.append(float(row["barrier"]))
            price.append(float(row["price"]))
            stderr.append(float(row["std_error"]))

    plt.figure()
    plt.plot(b, price, marker="o")
    plt.xlabel("Barrier B")
    plt.ylabel("MC price")
    plt.title("Barrier up-and-out call: price vs B")
    p1 = f"{out_prefix}_price.png"
    plt.savefig(p1, dpi=150, bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(b, stderr, marker="o", color="C1")
    plt.xlabel("Barrier B")
    plt.ylabel("MC std error")
    plt.title("Barrier sweep: statistical error vs B")
    p2 = f"{out_prefix}_stderr.png"
    plt.savefig(p2, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"wrote {p1}")
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
