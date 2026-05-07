#!/usr/bin/env python3
"""
HW4 Q4: reproduce-style Nusselt plots (Nu vs time / mesh) from convection.py CSVs.

  nu_history.csv format: header "time,Nu" then rows.

Usage (real data):
  python3 plot_q4_nusselt.py --figures-dir ../doc/figures \\
    --fig2 ../output/q4_ra1e2_n64/nu_history.csv \\
    --fig3 ../output/q4_ra1e2_n64/nu_history.csv:../output/q4_ra1e4_n64/nu_history.csv:... \\
    --fig3-labels 'Ra=10^2':'Ra=10^4':'Ra=10^5':'Ra=10^6' \\
    --fig4-final 16:path/nu.csv,32:path/nu.csv,64:path/nu.csv,128:path/nu.csv \\
    --fig5-ra1e4 16:path/nu.csv,32:path/nu.csv,64:path/nu.csv,128:path/nu.csv

Usage (demo curves qualitatively like the reference write-up; no Firedrake needed):
  python3 plot_q4_nusselt.py --figures-dir ../doc/figures --demo

Writes PNGs into figures-dir:
  hw4_q4_Nu_vs_t_Ra1e2_N64.png
  hw4_q4_Nu_vs_t_Ra_sweep_N64.png
  hw4_q4_Nu_final_vs_N_Ra1e4.png
  hw4_q4_Nu_vs_t_Ra1e4_mesh_sweep.png
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:  # pragma: no cover
    raise SystemExit("matplotlib is required: pip install matplotlib") from e


def load_nu_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    t_list: list[float] = []
    nu_list: list[float] = []
    with path.open() as f:
        # Accept both:
        # - convection.py output: CSV header "time,Nu"
        # - friend-style history: whitespace-separated with optional "# t  Nu" header
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if "," in s:
                parts = [p.strip() for p in s.split(",")]
            else:
                parts = s.split()
            if len(parts) < 2:
                continue
            try:
                t_list.append(float(parts[0]))
                nu_list.append(float(parts[1]))
            except ValueError:
                continue
    return np.asarray(t_list), np.asarray(nu_list)



def final_nu(t: np.ndarray, nu: np.ndarray, frac: float = 0.2) -> float:
    """Average Nu over the last `frac` of the time interval (robust to wiggles)."""
    if t.size == 0:
        return float("nan")
    t0 = float(t.max()) * (1.0 - frac)
    mask = t >= t0
    if not np.any(mask):
        return float(nu[-1])
    return float(np.nanmean(nu[mask]))


def style_axes(ax, xlabel: str, ylabel: str = r"$Nu$", *, blankenbach: bool = False) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.35)
    ax.axhline(
        1.0,
        color="k",
        lw=0.65,
        ls="--",
        alpha=0.45,
        label=r"$Nu=1$ (pure conduction)",
    )
    if blankenbach:
        ax.axhline(
            4.884,
            color="C2",
            lw=0.85,
            ls=":",
            alpha=0.75,
            label=r"Blankenbach $Nu\approx 4.884$",
        )


def synthetic_fig2() -> tuple[np.ndarray, np.ndarray]:
    t = np.linspace(0.0, 100_000.0, 4000)
    # Dip near t~15, recover to conductive Nu=1
    dip = 0.52 * np.exp(-0.5 * ((t - 18.0) / 14.0) ** 2)
    relax = (1.0 - np.exp(-t / 8000.0)) * 0.08 * np.sin(t / 900.0)
    nu = 1.0 - dip + relax
    return t, nu


def synthetic_fig3() -> list[tuple[str, np.ndarray, np.ndarray]]:
    curves = []
    t_full = np.linspace(0.0, 80_000.0, 3000)
    # Ra=10^2
    t2, n2 = synthetic_fig2()
    n2i = np.interp(t_full, t2, n2, left=n2[0], right=n2[-1])
    curves.append((r"$Ra=10^2$", t_full, n2i))
    # Ra=10^4 — oscillatory, mild overshoot
    rng = np.random.default_rng(0)
    nu4 = (
        1.0
        + 0.35 * np.sin(t_full / 2500.0) * (1.0 - np.exp(-t_full / 12_000.0))
        + 0.05 * rng.standard_normal(t_full.size)
    )
    curves.append((r"$Ra=10^4$", t_full, nu4))
    # Ra=10^5
    nu5 = 1.0 + 0.9 * np.sin(t_full / 1800.0 + 0.3) * (1.0 - np.exp(-t_full / 9000.0))
    curves.append((r"$Ra=10^5$", t_full, nu5))
    # Ra=10^6 — sharper transient, truncated (slow run)
    t6 = np.linspace(0.0, 35_000.0, 1200)
    nu6 = (
        1.0
        + 2.2 * np.sin(t6 / 800.0) * np.exp(-t6 / 25_000.0)
        - 0.8 * np.exp(-((t6 - 8000.0) / 2500.0) ** 2)
    )
    curves.append((r"$Ra=10^6$", t6, nu6))
    return curves


def synthetic_fig4() -> tuple[np.ndarray, np.ndarray]:
    Ns = np.array([16, 32, 64, 128], dtype=float)
    # Reference write-up: nearly flat near 1; add tiny drift for visibility
    nu = 1.0 + 0.02 * np.sin(Ns / 20.0)
    return Ns, nu


def synthetic_fig5() -> list[tuple[str, np.ndarray, np.ndarray]]:
    out: list[tuple[str, np.ndarray, np.ndarray]] = []
    base_t = np.linspace(0.0, 50_000.0, 2500)
    for i, N in enumerate([16, 32, 64, 128]):
        phase = 0.15 * i
        nu = (
            1.0
            + 0.25 * np.sin(base_t / 2200.0 + phase) * (1.0 - np.exp(-base_t / 10_000.0))
            + 0.03 * np.sin(base_t / 400.0 + N / 30.0)
        )
        out.append((rf"$N={N}$", base_t.copy(), nu))
    return out


def parse_kv_pairs(spec: str) -> list[tuple[int, Path]]:
    """'16:path/to.csv,32:path' -> [(16, Path), ...]"""
    out: list[tuple[int, Path]] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        k, v = chunk.split(":", 1)
        out.append((int(k.strip()), Path(v.strip())))
    return out


def _maybe_logx(ax, t: np.ndarray, nu: np.ndarray, logx: bool) -> tuple[np.ndarray, np.ndarray]:
    if not logx:
        return t, nu
    m = t > 0.0
    return t[m], nu[m]


def _robust_ylim(nu: np.ndarray) -> tuple[float, float] | None:
    """Pick y-limits that don't get dominated by transient spikes."""
    if nu.size < 5:
        return None
    nu = nu[np.isfinite(nu)]
    if nu.size < 5:
        return None
    lo = float(np.nanpercentile(nu, 1.0))
    hi = float(np.nanpercentile(nu, 99.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return None
    pad = 0.06 * (hi - lo)
    y0 = lo - pad
    y1 = hi + pad
    # Always include the conduction reference line.
    y0 = min(y0, 1.0)
    y1 = max(y1, 1.0)
    # Keep a non-degenerate window.
    if y1 - y0 < 1e-6:
        y0 -= 0.5
        y1 += 0.5
    return y0, y1


def plot_fig2(out: Path, t: np.ndarray, nu: np.ndarray, title: str, *, logx: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.25))
    t, nu = _maybe_logx(ax, t, nu, logx)
    ax.plot(t, nu, lw=1.4, color="C0", label=r"$Nu(t)$")
    ax.set_title(title)
    style_axes(ax, r"$t$", blankenbach=False)
    if logx:
        ax.set_xscale("log")
    yl = _robust_ylim(nu)
    if yl is not None:
        ax.set_ylim(*yl)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_fig3(
    out: Path, series: Iterable[tuple[str, np.ndarray, np.ndarray]], title: str, *, logx: bool = False
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    all_nu = []
    for i, (lab, t, nu) in enumerate(series):
        t, nu = _maybe_logx(ax, t, nu, logx)
        ax.plot(t, nu, lw=1.2, label=lab, color=f"C{i}")
        all_nu.append(nu)
    ax.set_title(title)
    style_axes(ax, r"$t$", blankenbach=True)
    if logx:
        ax.set_xscale("log")
    if all_nu:
        yl = _robust_ylim(np.concatenate(all_nu))
        if yl is not None:
            ax.set_ylim(*yl)
    h, l = ax.get_legend_handles_labels()
    if h:
        ax.legend(loc="best", fontsize=8, ncol=1)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_fig4(out: Path, Ns: np.ndarray, nu_fin: np.ndarray, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.25))
    ax.plot(Ns, nu_fin, "o-", lw=1.4, ms=7, color="C0")
    ax.axhline(4.884, color="C2", ls=":", lw=1.0, label=r"Blankenbach $4.884$")
    ax.set_xticks(Ns)
    ax.set_xticklabels([str(int(n)) for n in Ns])
    ax.set_title(title)
    ax.set_xlabel(r"Mesh size $N$ ($N\times N$ Q$_1$)")
    ax.set_ylabel(r"terminal / time-averaged $Nu$ (last 20% of $t$)")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_fig5(
    out: Path, series: Iterable[tuple[str, np.ndarray, np.ndarray]], title: str, *, logx: bool = False
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    all_nu = []
    for i, (lab, t, nu) in enumerate(series):
        t, nu = _maybe_logx(ax, t, nu, logx)
        ax.plot(t, nu, lw=1.1, label=lab, color=f"C{i}")
        all_nu.append(nu)
    ax.set_title(title)
    style_axes(ax, r"$t$", blankenbach=True)
    if logx:
        ax.set_xscale("log")
    if all_nu:
        yl = _robust_ylim(np.concatenate(all_nu))
        if yl is not None:
            ax.set_ylim(*yl)
    ax.legend(loc="best", fontsize=9, ncol=2)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-dir", type=Path, default=Path("../doc/figures"))
    ap.add_argument("--demo", action="store_true", help="Use synthetic curves (no CSVs); for layout / drafts.")
    ap.add_argument("--fig2", type=Path, help="CSV for Ra=10^2, N=64 Nu(t)")
    ap.add_argument("--fig3", type=str, help="Colon-separated CSV paths (order = labels)")
    ap.add_argument("--fig3-labels", type=str, default="", help="Colon-separated legend labels (optional)")
    ap.add_argument("--fig4-final", type=str, help="Comma list N:csvPath for Ra=10^4 terminal Nu vs N")
    ap.add_argument("--fig5-ra1e4", type=str, help="Comma list N:csvPath for Nu(t) at Ra=10^4")
    ap.add_argument("--logx", action="store_true", help="Log-scale the time axis for Nu(t) plots (drops t=0).")
    args = ap.parse_args()

    fig_dir: Path = args.figures_dir
    fig_dir.mkdir(parents=True, exist_ok=True)

    out2 = fig_dir / "hw4_q4_Nu_vs_t_Ra1e2_N64.png"
    out3 = fig_dir / "hw4_q4_Nu_vs_t_Ra_sweep_N64.png"
    out4 = fig_dir / "hw4_q4_Nu_final_vs_N_Ra1e4.png"
    out5 = fig_dir / "hw4_q4_Nu_vs_t_Ra1e4_mesh_sweep.png"

    if args.demo:
        t2, nu2 = synthetic_fig2()
        plot_fig2(out2, t2, nu2, r"$Nu(t)$ for $Ra=10^2$, $64\times 64$ (demo curve)", logx=args.logx)
        plot_fig3(
            out3,
            [(lab, tt, nn) for lab, tt, nn in synthetic_fig3()],
            r"$Nu(t)$ for $Ra\in\{10^2,10^4,10^5,10^6\}$, $64\times 64$ (demo)",
            logx=args.logx,
        )
        N4, nu4 = synthetic_fig4()
        plot_fig4(
            out4,
            N4,
            nu4,
            r"Terminal $Nu$ vs mesh for $Ra=10^4$ (demo; replace with simulation CSVs)",
        )
        plot_fig5(
            out5,
            [(lab, tt, nn) for lab, tt, nn in synthetic_fig5()],
            r"$Nu(t)$ for $Ra=10^4$ and mesh sizes $N\in\{16,32,64,128\}$ (demo)",
            logx=args.logx,
        )
        print(f"Wrote demo PNGs to {fig_dir.resolve()}")
        return

    if args.fig2:
        t, nu = load_nu_csv(args.fig2)
        plot_fig2(out2, t, nu, r"$Nu(t)$ for $Ra=10^2$, $64\times 64$", logx=args.logx)
    if args.fig3:
        paths = [Path(p) for p in args.fig3.split(":") if p.strip()]
        labels = [p.strip() for p in args.fig3_labels.split(":")] if args.fig3_labels else []
        if len(labels) != len(paths):
            labels = [Path(p).parent.name for p in paths]
        series = []
        for lab, p in zip(labels, paths):
            t, nu = load_nu_csv(p)
            series.append((lab, t, nu))
        plot_fig3(out3, series, r"$Nu(t)$ Rayleigh sweep, $64\times 64$", logx=args.logx)
    if args.fig4_final:
        pairs = parse_kv_pairs(args.fig4_final)
        pairs.sort(key=lambda x: x[0])
        Ns = np.array([n for n, _ in pairs], dtype=float)
        nfs = []
        for n, p in pairs:
            t, nu = load_nu_csv(p)
            nfs.append(final_nu(t, nu))
        plot_fig4(
            out4,
            Ns,
            np.asarray(nfs),
            r"Terminal $Nu$ vs mesh for $Ra=10^4$",
        )
    if args.fig5_ra1e4:
        pairs = parse_kv_pairs(args.fig5_ra1e4)
        pairs.sort(key=lambda x: x[0])
        series = []
        for n, p in pairs:
            t, nu = load_nu_csv(p)
            series.append((rf"$N={n}$", t, nu))
        plot_fig5(out5, series, r"$Nu(t)$ for $Ra=10^4$ (mesh refinement)", logx=args.logx)

    print(f"Done. Outputs in {fig_dir.resolve()}")


if __name__ == "__main__":
    main()
