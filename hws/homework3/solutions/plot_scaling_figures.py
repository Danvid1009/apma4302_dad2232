#!/usr/bin/env python3
"""
Two-panel scaling figure for HW3 Problem 8 (matches course-style plots).
Reads solutions/logs/scaling_hw3_runs.txt if present; else uses embedded sample data.
Output: solutions/figures/scaling_runtime_error.png (and .pdf)
"""
from __future__ import annotations

import os
import re
import sys

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Need matplotlib and numpy: pip install matplotlib numpy", file=sys.stderr)
    sys.exit(1)

LOG = os.path.join(os.path.dirname(__file__), "logs", "scaling_hw3_runs.txt")
OUT_DIR = os.path.join(os.path.dirname(__file__), "figures")


def parse_log(path: str) -> tuple[dict[int, list[float]], list[int], list[float]]:
    """Return times[np] lists aligned with Ns, plus Ns and rel_err per grid."""
    times: dict[int, list[float]] = {1: [], 2: [], 4: []}
    Ns: list[int] = []
    errs: list[float] = []
    if not os.path.isfile(path):
        return times, Ns, errs
    refine_N = {2: 33, 3: 65, 4: 129, 5: 257, 6: 513}
    by_refine: dict[int, dict] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            m = re.match(
                r"refine=(\d+)\s+np=(\d+)\s+real=([\d.]+)\s+SNES=\d+\s+rel_err=([\d.eE+-]+)",
                line,
            )
            if not m:
                continue
            r, np_, t, e = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
            by_refine.setdefault(r, {"t": {}, "e": e})
            by_refine[r]["t"][np_] = t
            by_refine[r]["e"] = e
    for r in sorted(by_refine):
        if r not in refine_N:
            continue
        td = by_refine[r]["t"]
        if not all(k in td for k in (1, 2, 4)):
            continue
        Ns.append(refine_N[r])
        errs.append(by_refine[r]["e"])
        for np_ in (1, 2, 4):
            times[np_].append(td[np_])
    return times, Ns, errs


def default_data():
    """Embedded data from scaling_hw3_runs.txt (sample run)."""
    Ns = [33, 65, 129, 257, 513]
    errs = [8.09754e-4, 2.02508e-4, 5.06487e-5, 1.26657e-5, 3.16692e-6]
    times = {
        1: [0.25, 0.19, 0.25, 0.52, 1.55],
        2: [0.21, 0.19, 0.27, 0.50, 1.43],
        4: [0.21, 0.22, 0.33, 0.51, 1.39],
    }
    return times, Ns, errs


def main():
    times, Ns, errs = parse_log(LOG)
    if not Ns:
        times, Ns, errs = default_data()

    N = np.array(Ns, dtype=float)
    err = np.array(errs, dtype=float)

    # O(h^2) ~ 1/N^2 reference anchored at finest grid
    N0, e0 = N[-1], err[-1]
    err_ref = e0 * (N0 / N) ** 2

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    colors = {1: "C0", 2: "C1", 4: "C2"}
    for np_ in (1, 2, 4):
        y = np.array(times[np_], dtype=float)
        ax1.plot(N, y, "o-", color=colors[np_], lw=1.5, ms=7, label=f"{np_} proc(s)")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Grid size ($N \times N$)")
    ax1.set_ylabel("Time (sec)")
    ax1.set_title("Run time vs grid size")
    ax1.grid(True, which="both", ls="--", alpha=0.35)
    ax1.legend(loc="upper left")
    ax1.set_xticks(N)
    ax1.set_xticklabels([f"{int(n)}$\\times${int(n)}" for n in N], rotation=35, ha="right")

    ax2.plot(N, err, "o-", color="C0", lw=1.5, ms=7, label="computed")
    ax2.plot(N, err_ref, "k--", lw=1.5, label=r"$O(h^2)$ reference ($\propto N^{-2}$)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"Grid size ($N \times N$)")
    ax2.set_ylabel("Relative error")
    ax2.set_title("Relative error vs grid size")
    ax2.grid(True, which="both", ls="--", alpha=0.35)
    ax2.legend(loc="upper right")
    ax2.set_xticks(N)
    ax2.set_xticklabels([f"{int(n)}$\\times${int(n)}" for n in N], rotation=35, ha="right")

    os.makedirs(OUT_DIR, exist_ok=True)
    png = os.path.join(OUT_DIR, "scaling_runtime_error.png")
    pdf = os.path.join(OUT_DIR, "scaling_runtime_error.pdf")
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    print("Wrote", png, "and", pdf)


if __name__ == "__main__":
    main()
