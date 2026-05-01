#!/usr/bin/env python3
"""Extra credit (HW4 item 1): mesh sweep for Ra in {1e4,1e5,1e6} × N grid, using convection.py (DAE).

Writes summary CSV for tables / plots. Run inside Firedrake+firedrake-ts (same as Q4)."""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def last_nu(csv_path: Path) -> float:
    with csv_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return float("nan")
    return float(rows[-1]["Nu"])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--t-max", type=float, default=5000.0, help="Per-case end time (raise for steadier Nu)")
    p.add_argument("--dt", type=float, default=0.1)
    p.add_argument(
        "--nu-every",
        type=int,
        default=0,
        help="Passed to convection.py (0 = use script default); larger = fewer Nu logs, faster",
    )
    p.add_argument(
        "--output-root",
        type=Path,
        default=Path("output/bonus_convergence"),
        help="Directory under cwd (usually hw4/python)",
    )
    p.add_argument(
        "--ra",
        type=float,
        nargs="+",
        default=[1e4, 1e5, 1e6],
    )
    p.add_argument("--n", type=int, nargs="+", default=[16, 32, 64, 128])
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    root: Path = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).resolve().parent / "convection.py"
    py = sys.executable

    summary_path = root / "summary_nu_final.csv"
    summary_rows: list[tuple[float, int, float, str]] = []

    for ra in args.ra:
        for n in args.n:
            case = root / f"Ra{ra:g}_N{n}"
            case.mkdir(parents=True, exist_ok=True)
            cmd = [
                py,
                str(script),
                "--ra",
                str(ra),
                "--n",
                str(n),
                "--t-max",
                str(args.t_max),
                "--dt",
                str(args.dt),
                "--output-dir",
                str(case),
                "--vtk-every",
                "0",
            ]
            if args.nu_every > 0:
                cmd.extend(["--nu-every", str(args.nu_every)])
            print("+", " ".join(cmd))
            if args.dry_run:
                continue
            subprocess.run(cmd, check=True)
            nu_path = case / "nu_history.csv"
            nu_f = last_nu(nu_path) if nu_path.exists() else float("nan")
            summary_rows.append((ra, n, nu_f, str(nu_path)))
            print(f"  -> Nu_final ≈ {nu_f} ({nu_path})")

    if not args.dry_run:
        with summary_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Ra", "N", "Nu_final", "case_dir_nu_history"])
            w.writerows(summary_rows)
        print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
