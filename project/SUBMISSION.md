# APMA 4302 final project — submission bundle

## What to submit

| Artifact | Path |
| --- | --- |
| **Report (PDF)** | [`submission/4302_Final_Project_Submission.pdf`](submission/4302_Final_Project_Submission.pdf) |
| **LaTeX source** | [`overleaf/main.tex`](overleaf/main.tex) + [`overleaf/figures/`](overleaf/figures/) |
| **Code & scripts** | [`src/`](src/), [`scripts/`](scripts/) |
| **Repro notes** | [`README.md`](README.md), [`PROJECT_OUTLINE.md`](PROJECT_OUTLINE.md), [`docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md`](docs/PROJECT_SUMMARY_AND_PETSC_CHECKLIST.md) |

**Dropped locally (regenerate as needed):** `project/.venv/`, all `output/paraview/vtp/*.vtp` (~12 MB+), duplicate root PDF, `overleaf.zip`. Only `submission/4302_Final_Project_Submission.pdf` is kept as the canonical report.

## Rebuild the PDF

```bash
cd overleaf
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

See [`overleaf/README_COMPILE.md`](overleaf/README_COMPILE.md).

## Branch

Course submission snapshot: git branch **`final-project`**.
