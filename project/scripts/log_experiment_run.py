"""
Append one row to output/registry/experiment_runs.csv for bookkeeping.

Use after scaling, convergence, matrix batches, or cluster jobs so outputs stay traceable.
"""

import argparse
import csv
import datetime as dt
import pathlib
import socket
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Append a row to the HPC experiment run registry CSV.")
    p.add_argument("--key", required=True, help="Short experiment id, e.g. strong_european_p8.")
    p.add_argument("--ranks", default="", help="MPI rank list or description, e.g. 1;2;4 or 8.")
    p.add_argument(
        "--output",
        action="append",
        dest="outputs",
        default=[],
        help="Output artifact (repeatable): CSV, PNG, md path, etc.",
    )
    p.add_argument("--status", default="completed", help="completed | failed | archived | running")
    p.add_argument("--notes", default="", help="Free text; keep commas minimal or quote manually in CSV later.")
    p.add_argument("--cmd", default="", help="Optional shell command string for reproducibility.")
    p.add_argument(
        "--registry-dir",
        type=str,
        default="output/registry",
        help="Directory containing experiment_runs.csv",
    )
    return p.parse_args()


def main():
    args = parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    reg_dir = (root / args.registry_dir).resolve()
    reg_dir.mkdir(parents=True, exist_ok=True)
    csv_path = reg_dir / "experiment_runs.csv"

    fieldnames = [
        "timestamp_utc",
        "experiment_key",
        "hostname",
        "ranks",
        "output_paths",
        "status",
        "notes",
        "command",
    ]
    row = {
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "experiment_key": args.key,
        "hostname": socket.gethostname(),
        "ranks": args.ranks,
        "output_paths": ";".join(args.outputs) if args.outputs else "",
        "status": args.status,
        "notes": args.notes.replace("\n", " ").strip(),
        "command": args.cmd.replace("\n", " ").strip(),
    }

    new_file = not csv_path.exists()
    with csv_path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if new_file:
            writer.writeheader()
        writer.writerow(row)

    print(f"appended -> {csv_path.relative_to(root)}")
    print("tip: mirror a one-line summary in docs/HPC_EXPERIMENT_REGISTRY.md if it is a milestone run.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
