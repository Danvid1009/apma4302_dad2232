# Homework 4 write-up

**Main submission PDF:** `../HW4_4302.pdf`  

The solutions write-up is contained in that PDF. Everything else (source LaTeX, code, figures, and raw Q4 histories used for plots) is in the surrounding `homework4/` directory.

---

# Building the PDF from LaTeX

From the **`homework4/doc`** directory (so relative paths to `figures/` resolve):

```bash
cd doc
pdflatex -interaction=nonstopmode solutions.tex
pdflatex -interaction=nonstopmode solutions.tex   # second pass for references
```

- **`solutions.tex`** — main submission write-up (includes handout pages from `figures/hw4_*.png`).
- **`hw4.tex`** — original assignment text only (optional compile).

Outputs:

- `solutions.pdf` (rename to `hw4.pdf` if you want a single hand-in name: `cp solutions.pdf hw4.pdf`).

The committed **`../HW4_4302.pdf`** is a recent build of the solutions document; rebuild after you update text or regenerate figures under `figures/`.
