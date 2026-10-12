# CEB-ViTAL manuscript (LaTeX)

The authors are Ali Hamza and Javed Iqbal. The document uses the IEEEtran journal format.

## VS Code (LaTeX Workshop)

1. Reload the window once so the extension activates.
2. Open `main.tex` and save it. The recipe in `../.vscode/settings.json` builds the PDF and opens it in a tab. The recipe runs pdflatex, then bibtex, then pdflatex twice.
3. You can also build manually:

```bash
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Overleaf

On Overleaf, choose New Project, then Upload Project, and pick `../CEB-ViTAL_Overleaf.zip`. Set the main document to `main.tex`; Overleaf compiles it as is.

## Where the numbers come from

Do not edit `numbers.tex`, `numbers_extra.tex`, `tables/*.tex` or `figures/*.pdf` by hand. They are generated from the simulation data:

```bash
cd ../paper && python3 make_latex_assets.py && python3 make_paper_figures.py && python3 make_extra_figures.py
```

## Before submission

Search the PDF for red `[...]` placeholders and fill each one:

- affiliations and e-mail addresses
- corresponding author
- journal name
- repository URL
- AI-use declaration
- author contributions

Also extend the related-work section with your own literature search, and check every `refs.bib` entry against the publisher's record.
