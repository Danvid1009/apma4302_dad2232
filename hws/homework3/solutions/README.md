# APMA 4302 — Homework 3 (solutions workspace)

Use this folder for write-ups, plots, option files, and notes. Implement the code in **`p4pdes`** (`c/ch4/reaction.c` → executable **`reaction2d`**); reference **`apma4302-methods/homework/hw3/poisson2d/`** for the linear Poisson template and VTK patterns.

---

## Part 1 (5 pts) — Discrete nonlinear system

**Deliverable:** Written derivation.

- Discretize \(-\nabla^2 u + \gamma u^p = f\) on a 2D finite-difference grid (unit square).
- Write the **nonlinear algebraic system** in the form **\(F(\mathbf{u}) = 0\)** (stack interior + boundary unknowns as required by your convention).
- Include **auxiliary equations** (or explicit treatment) for **non-homogeneous Dirichlet** data: \(u = u_{\text{exact}}\) on the boundary.
- State clearly how indices map to grid points and how the Laplacian stencil enters \(F\).

---

## Part 2 (5 pts) — Jacobian

**Deliverable:** Symbolic formulas for **\(J_{i,j}(\mathbf{u}) = \partial F_i / \partial u_j\)**.

- Differentiate your discrete \(F\) with respect to each unknown \(u_j\).
- Account for boundary equations (rows where \(u\) is fixed) so \(J\) is well-defined for Newton.
- You may use a sparse structure description (which entries are nonzero, 5-point + diagonal from \(\gamma p u^{p-1}\)).

---

## Part 3 (20 pts) — `reaction2d` with PETSc SNES

**Deliverable:** Working parallel code + VTK output.

### 3(a) CLI

Run like:

```bash
mpirun -n <num_procs> ./reaction2d -rct_p <p> -rct_gamma <gamma> [-rct_linear_f]
```

plus usual PETSc options (`-da_refine`, KSP/SNES, etc.).

- **`-rct_linear_f`** (boolean): if **true**, use  
  \(f(x,y) = -\nabla^2 u_{\text{exact}}\)  
  if **false**, use the full MMS RHS  
  \(f(x,y) = -\nabla^2 u_{\text{exact}} + \gamma u_{\text{exact}}^p\).
- Reuse / adapt **`poisson2d.c`** logic for \(u_{\text{exact}}\) and its derivatives.

### 3(b) Boundary conditions

- On **all** boundaries: \(u = u_{\text{exact}}\), **regardless** of `-rct_linear_f`.

### 3(c) Residual and Jacobian

- Implement **\(F(\mathbf{u})\)** and a user-assembled **\(J(\mathbf{u})\)** (or `MatShell` if allowed by your design—handout expects explicit Jacobian capability).

### 3(d) Finite-difference Jacobians

- Code must run correctly with PETSc’s **matrix-free / FD** Jacobian options, e.g.  
  `-snes_fd`, `-snes_fd_color`, `-snes_mf`, `-snes_mf_operator`.
- Verify with **`-snes_test_jacobians`** against your hand-coded \(J\).

### 3(e) VTK

- Write the numerical solution in **VTK** format for **ParaView** (consistent with the style used in `poisson2d.c`).

---

## Part 4 (5 pts) — Linear limit \(\gamma = 0\), \(65\times 65\)

**Deliverable:** Short report + **`options_file_gamma0`**.

- Grid: **\(65 \times 65\)** (choose the `-da_refine` / sizing options that produce this).
- Show **nonlinear SNES** and **linear** formulations agree to expected discretization / tolerance error when \(\gamma = 0\).
- Tune **linear solver** and tolerances so SNES reaches **\(\|F\|_2 < 10^{-10}\)** in **one** Newton iteration.
- Report:
  - **KSP iterations** at each Newton step (with `-ksp_monitor`),
  - final **relative error** \(\|\mathbf{u} - \mathbf{u}_{\text{exact}}\|_2 / \|\mathbf{u}_{\text{exact}}\|_2\).
- Hints: `-snes_monitor`, `-ksp_rtol`, `-ksp_atol`; try preconditioners / direct solvers; must run **in parallel**.
- Save your chosen PETSc options in **`options_file_gamma0`** and submit it.

---

## Part 5 (2 pts) — Finite-difference Jacobian comparison

**Deliverable:** Notes + **`options_file_fd`**.

- Run with a **finite-difference Jacobian**; confirm you obtain the **same solution** (within tolerance) as the analytic \(J\).
- Describe any differences in **convergence** or **final error**.
- Put the exact command-line (or options file contents) in **`options_file_fd`**.
- Compare **runtimes** (user Jacobian vs FD Jacobian) using **`-log_view`**.

---

## Part 6 (5 pts) — Nonlinear case \(\gamma = 100\), \(p = 3\)

**Deliverable:** Report + residual plot.

- Use the **same SNES/KSP stopping criteria** as in Part 4 (unless you document a deliberate change).
- Report:
  - total **Newton** iterations,
  - **KSP iterations per Newton** step,
  - final **relative error** \(\|\mathbf{u} - \mathbf{u}_{\text{exact}}\|_2 / \|\mathbf{u}_{\text{exact}}\|_2\).
- **Plot** \(\|F(\mathbf{u})\|_2\) vs Newton iteration (convergence history).

---

## Part 7 (5 pts) — ParaView, `-rct_linear_f`

**Deliverable:** Screenshot + short discussion.

- Rerun with **`-rct_linear_f`** for \(\gamma = 100\), \(p = 3\).
- Visualize in **ParaView**; include a **screenshot** in your submission.
- **Describe** how the nonlinearity changes the solution compared to the **linear** case \(\gamma = 0\).

---

## Part 8 (10 pts) — Scaling study

**Deliverable:** Table and/or plots + commentary.

For the **nonlinear** problem, compare (use **`-log_view`** for times):

| Quantity | What to record |
|----------|----------------|
| Run time | Wall time / PETSc event log |
| Newton iters | Count |
| Relative error | \(\|\mathbf{u}-\mathbf{u}_{\text{exact}}\|_2/\|\mathbf{u}_{\text{exact}}\|_2\) |

**Parameter sweep:**

- **`-da_refine`** \(= 2,3,4,5,6\) → grids **\(33^2\)** through **\(513^2\)** (as in the handout).
- **MPI ranks:** **1, 2, 4**.

Comment on scaling in **problem size** and in **processor count** (strong scaling trends, expected limitations).

---

## Suggested layout inside `solutions/`

| Path | Purpose |
|------|---------|
| `writeup/` | LaTeX/PDF for Parts 1–2, narrative for 4–8 |
| `figures/` | Plots (residual history, scaling, ParaView exports) |
| `options_file_gamma0` | PETSc options for Part 4 |
| `options_file_fd` | PETSc options for FD Jacobian runs |
| `logs/` | Saved `-log_view` text for timing tables |

Adjust names to match your course submission rules.
