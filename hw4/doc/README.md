# Building `hw4.pdf` from LaTeX

From the **`hw4/doc`** directory (so relative paths to `figures/` resolve):

```bash
cd hw4/doc
pdflatex -interaction=nonstopmode solutions.tex
pdflatex -interaction=nonstopmode solutions.tex   # second pass for references
```

- **`solutions.tex`** — main submission write-up (includes handout pages from `figures/hw4_*.png`).
- **`hw4.tex`** — original assignment text only (optional compile).

Outputs:

- `solutions.pdf` (rename to `hw4.pdf` if you want a single hand-in name: `cp solutions.pdf hw4.pdf`).

The committed **`hw4.pdf`** is a recent build of the solutions document; rebuild after you replace placeholders or add figures under `figures/`.
