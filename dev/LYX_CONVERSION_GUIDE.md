# Converting Research.tex to Research.lyx

## Quick Method

1. Open LyX 2.4+
2. File → Import → LaTeX (plain)
3. Select `Research.tex`
4. LyX will parse the LaTeX and generate the `.lyx` file

## Post-Import Adjustments

After importing, you may need to:

1. **Document Class**: Verify LyX recognized `article` class (Document → Settings → Document Class)
2. **Theorem Environments**: The `\newtheorem` definitions will be imported as ERT (raw LaTeX). In LyX, go to Document → Settings → Modules and add:
   - "AMS Theorems" module (provides Theorem, Proposition, Lemma, Corollary, Remark)
   - Then convert ERT theorem blocks to native LyX theorem insets
3. **Hypothesis Environment**: The custom `\newtheorem{hyp}{Hypothesis}` with `\renewcommand{\thehyp}{V\arabic{hyp}}` will remain as ERT. This is fine — LyX handles custom environments via ERT gracefully.
4. **Bibliography**: Ensure `Research.bib` is in the same directory. Insert → List/TOC → BibTeX Bibliography and point to `Research.bib`.
5. **Figure Paths**: The tex uses relative paths `figures/figure_N_name`. In LyX, re-insert figures via Insert → Graphics and browse to the files, or verify the paths resolve correctly.

## File Organization

Place these files in your `report/` directory:
```
report/
├── Research.tex          (revised source)
├── Research.lyx          (after LyX import)
├── Research.bib          (bibliography, unchanged)
└── figures/              (symlink or copy from ../figures/output/)
    ├── figure_1_admissibility_binding.pdf
    ├── figure_2_peak_configurations.pdf
    ├── figure_3_alpha_vs_theta_vs_mu.pdf
    ├── figure_4_profile_decomposition.pdf
    ├── figure_5_concentration_vs_peaks.pdf
    ├── figure_6_interaction_decay.pdf
    └── figure_7_methodology_flowchart.pdf
```

## Compiling the TeX Directly

```bash
cd report
pdflatex Research.tex
bibtex Research
pdflatex Research.tex
pdflatex Research.tex
```
