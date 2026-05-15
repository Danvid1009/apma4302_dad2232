# APMA 4302 final project — submission bundle

## What to submit

| Artifact | Path |
| --- | --- |
| **Report (PDF)** | [`submission/4302_Final_Project.pdf`](submission/4302_Final_Project.pdf) |
| **LaTeX source** | [`overleaf/main.tex`](overleaf/main.tex) + [`overleaf/figures/`](overleaf/figures/) |
| **Code & scripts** | [`src/`](src/), [`scripts/`](scripts/) |
| **Repro notes** | [`README.md`](README.md), [`PROJECT_OUTLINE.md`](PROJECT_OUTLINE.md), [`docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md`](docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md) |

The PDF at the project root (`4302_Final_Project.pdf`) is a copy of the same file in `submission/` for convenience.

## Rebuild the PDF

```bash
cd overleaf
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

See [`overleaf/README_COMPILE.md`](overleaf/README_COMPILE.md).

## Branch

Course submission snapshot: git branch **`final-project`**.
