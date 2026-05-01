import argparse
import csv
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from montecarlo.gbm import simulate_gbm_paths


def parse_args():
    parser = argparse.ArgumentParser(description="Export Monte Carlo paths for ParaView.")
    parser.add_argument("--s0", type=float, default=100.0)
    parser.add_argument("--r", type=float, default=0.05)
    parser.add_argument("--sigma", type=float, default=0.2)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--paths", type=int, default=200)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--barrier", type=float, default=130.0)
    parser.add_argument("--out", type=str, default="output/path_bundle.csv")
    parser.add_argument("--barrier-out", type=str, default="output/barrier_plane.csv")
    return parser.parse_args()


def export_paths(out_path: pathlib.Path, paths: np.ndarray, t: float):
    n_paths, n_cols = paths.shape
    times = np.linspace(0.0, t, n_cols)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["path_id", "time", "stock_price"])
        for i in range(n_paths):
            for j, tt in enumerate(times):
                writer.writerow([i, float(tt), float(paths[i, j])])


def export_barrier_plane(out_path: pathlib.Path, barrier: float, t: float, n_paths: int):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["path_id", "time", "barrier"])
        for i in range(n_paths):
            writer.writerow([i, 0.0, barrier])
            writer.writerow([i, t, barrier])


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    paths = simulate_gbm_paths(args.s0, args.r, args.sigma, args.t, args.steps, args.paths, rng)

    export_paths(pathlib.Path(args.out), paths, args.t)
    export_barrier_plane(pathlib.Path(args.barrier_out), args.barrier, args.t, args.paths)

    print(f"wrote paths to {args.out}")
    print(f"wrote barrier plane helper to {args.barrier_out}")


if __name__ == "__main__":
    main()
