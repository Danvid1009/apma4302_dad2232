"""
Batch strong/weak scaling runs across option types (and optional antithetic flag).

Writes one CSV per (mode, option[, antithetic]) plus a small JSON manifest for
downstream plotting and report tables.
"""

import argparse
import json
import pathlib
import subprocess
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser(description="Run a matrix of scaling experiments.")
    parser.add_argument("--modes", nargs="+", choices=["strong", "weak"], default=["strong", "weak"])
    parser.add_argument(
        "--options",
        nargs="+",
        choices=["european", "asian", "barrier"],
        default=["european", "asian", "barrier"],
    )
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--base-paths", type=int, default=200000)
    parser.add_argument("--steps", type=int, default=252)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--antithetic", action="store_true", help="Also run an antithetic matrix (adds runs).")
    parser.add_argument("--out-dir", type=str, default="output/experiments")
    return parser.parse_args()


def run_scaling(mode: str, option: str, antithetic: bool, args):
    tag = "anti" if antithetic else "plain"
    out = pathlib.Path(args.out_dir) / f"scaling_{mode}_{option}_{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "scripts/run_scaling.py",
        "--mode",
        mode,
        "--option",
        option,
        "--ranks",
        *[str(p) for p in args.ranks],
        "--base-paths",
        str(args.base_paths),
        "--steps",
        str(args.steps),
        "--seed",
        str(args.seed),
        "--out",
        str(out),
    ]
    if antithetic:
        cmd.append("--antithetic")
    t0 = time.perf_counter()
    subprocess.run(cmd, check=True, cwd=pathlib.Path(__file__).resolve().parents[1])
    elapsed = time.perf_counter() - t0
    return out, elapsed


def main():
    args = parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    manifest = {
        "started_unix_s": time.time(),
        "runs": [],
    }
    anti_flags = [False, True] if args.antithetic else [False]
    for mode in args.modes:
        for option in args.options:
            for anti in anti_flags:
                rel_csv, elapsed = run_scaling(mode, option, anti, args)
                manifest["runs"].append(
                    {
                        "mode": mode,
                        "option": option,
                        "antithetic": anti,
                        "csv": str(rel_csv.resolve().relative_to(root.resolve())),
                        "wall_s": elapsed,
                    }
                )
                print(f"done {mode} {option} anti={anti} -> {rel_csv} ({elapsed:.1f}s)")

    manifest["finished_unix_s"] = time.time()
    man_path = pathlib.Path(args.out_dir) / "experiment_matrix_manifest.json"
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {man_path}")
    print(
        "tip: log this batch to the registry, e.g.\n"
        f"  python scripts/log_experiment_run.py --key experiment_matrix "
        f'--notes "see {man_path}" --output {man_path}'
    )


if __name__ == "__main__":
    main()
