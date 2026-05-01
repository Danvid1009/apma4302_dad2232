#!/usr/bin/env python3
"""Plot ||F|| vs SNES iteration from Part 6 (edit norms list if you re-run with -snes_monitor)."""
import os
import sys

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Install matplotlib: pip install matplotlib", file=sys.stderr)
    sys.exit(1)

# Values from: mpirun ... ./reaction2d -rct_gamma 100 -rct_p 3 -da_refine 3 \
#   -snes_atol 1e-12 -snes_rtol 0 -ksp_type preonly -pc_type lu -snes_monitor
SNES_FNORMS = [
    2.164823134561e-04,
    1.898785195959e-08,
    7.203864841333e-15,
]
it = list(range(len(SNES_FNORMS)))

out_dir = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(out_dir, exist_ok=True)
out_png = os.path.join(out_dir, "snes_residual_history.png")

fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(it, SNES_FNORMS, "o-", lw=1.5, ms=8)
ax.set_xlabel("SNES iteration")
ax.set_ylabel(r"$\|F(\mathbf{u})\|_2$ (printed SNES function norm)")
ax.set_title(r"Nonlinear residual history ($\gamma=100$, $p=3$, $65\times 65$)")
ax.grid(True, which="both", ls="--", alpha=0.4)
fig.tight_layout()
fig.savefig(out_png, dpi=150)
print("Wrote", out_png)
