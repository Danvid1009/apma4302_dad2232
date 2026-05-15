# LaTeX report — compile instructions

## Where the sources live

| Path | Role |
|------|------|
| `apma4302_dad2232/project/overleaf/main.tex` | **Main document** (article class, `thebibliography` — no BibTeX required) |
| `apma4302_dad2232/project/overleaf/figures/` | **PNG figures** synced from `../visuals/` and `../output/` for a self-contained zip |

## Local compile (Terminal)

```bash
cd /Users/dan/Desktop/Columbia/HPC_4302/apma4302_dad2232/project/overleaf
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Output: `main.pdf` in the same folder. Run **twice** so the table of contents and cross-references settle.

### Refresh figures before compiling

From the **project root** (`.../project/`, parent of `overleaf/`):

```bash
mkdir -p overleaf/figures
cp -f visuals/*.png overleaf/figures/   # broad refresh; or cherry-pick as in repo scripts
cp -f output/convergence.png output/error_surface_abs_err.png overleaf/figures/ 2>/dev/null || true
```

(Use a curated `cp` list if the glob is too large.)

## Overleaf (browser)

1. Zip the **`overleaf/`** folder **including** `figures/` and upload **Upload Project**.
2. Set the Main document to **`main.tex`** (Menu → Main file).
3. Recompiler: **pdfLaTeX** (no BibTeX step unless you later switch to `.bib`).
4. If uploads omit `figures/`, use the **project menu → Add files** and add PNGs under `figures/` with **exact filenames** referenced in `main.tex`.

## Remaining submission steps (final polish)

See **`main.tex`**, subsection ``Next steps: low-hanging fruit'' (LaTeX label `sec:remaining`): keep `figures/` synced whenever you regenerate plots, update the author line, check registry/numbers consistency, optional quotes comparison, and run `pdflatex` twice. The main PDF already embeds representative ParaView / trajectory stills; refresh those PNGs if you change export parameters.

To add **more** ParaView panels later: new PNGs under `figures/` + new `figure` blocks in `main.tex`; workflow: **`PROJECT_OUTLINE.md` Appendix A** (repo root) or stub `PARAVIEW_GUIDE_PETSC.md`.

## Primary references cited in `main.tex`

- **E. Bueler**, *PETSc for Partial Differential Equations: Numerical Solutions in C and Python*, SIAM, **2020** — main PDE / PETSc numerics text.
- **PETSc** user manual + original PETSc paper (Balay et al.).
- **Black–Scholes / Merton** for options.
- **ParaView** chapter (Ahrens et al.).
