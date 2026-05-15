import argparse
import csv
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]


def parse_args():
    parser = argparse.ArgumentParser(description="Run PETSc/MPI strong/weak scaling experiments.")
    parser.add_argument("--mode", choices=["strong", "weak"], required=True)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--base-paths", type=int, default=200000)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--option", choices=["european", "asian", "barrier"], default="european")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--antithetic",
        action="store_true",
        help="Forward --antithetic to run_pricing.py (paths forced even).",
    )
    parser.add_argument("--out", type=str, default="output/scaling_results.csv")
    return parser.parse_args()


def run_once(rank_count: int, option: str, paths: int, steps: int, seed: int, antithetic: bool) -> float:
    cmd = [
        "mpirun",
        "-n",
        str(rank_count),
        sys.executable,
        str(ROOT / "scripts/run_pricing.py"),
        "--option",
        option,
        "--paths",
        str(paths),
        "--steps",
        str(steps),
        "--seed",
        str(seed),
    ]
    if antithetic:
        cmd.append("--antithetic")
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, cwd=ROOT)
    t1 = time.perf_counter()
    return t1 - t0


def main():
    args = parse_args()
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    t1_baseline = None

    for p in args.ranks:
        if args.mode == "strong":
            total_paths = args.base_paths
            paths_per_rank = max(1, total_paths // p)
        else:
            paths_per_rank = args.base_paths
            total_paths = paths_per_rank * p

        elapsed = run_once(
            rank_count=p,
            option=args.option,
            paths=total_paths,
            steps=args.steps,
            seed=args.seed,
            antithetic=args.antithetic,
        )

        if p == 1:
            t1_baseline = elapsed
        speedup = (t1_baseline / elapsed) if t1_baseline else 1.0
        efficiency = speedup / p

        rows.append(
            {
                "mode": args.mode,
                "ranks": p,
                "total_paths": total_paths,
                "paths_per_rank": paths_per_rank,
                "elapsed_s": elapsed,
                "speedup": speedup,
                "efficiency": efficiency,
            }
        )
        print(
            f"mode={args.mode} p={p} total_paths={total_paths} elapsed={elapsed:.4f}s "
            f"speedup={speedup:.3f} efficiency={efficiency:.3f}"
        )

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "mode",
                "ranks",
                "total_paths",
                "paths_per_rank",
                "elapsed_s",
                "speedup",
                "efficiency",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
