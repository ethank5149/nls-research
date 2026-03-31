# Comprehensive Improvement Report
## Taking "Existence, Uniqueness, and Stability of Ground-State Solutions for the NLS in Optical Fibers" to Publication Quality

**Prepared for:** Ethan Knox & Eduard Kirr  
**Date:** March 30, 2026  
**Scope:** Full audit of the current `Research.tex` report, repository structure, mathematical content, figures, and presentation, with prioritized recommendations for achieving publication in a peer-reviewed mathematics or mathematical physics journal (e.g., *Communications in Mathematical Physics*, *Journal of Functional Analysis*, *Archive for Rational Mechanics and Analysis*, *Nonlinearity*).

---

## Executive Summary

The current report presents genuinely original mathematical work: a systematic classification of admissible multi-peak configurations for the stationary NLS with external potential, and a carefully structured attack on the kernel obstruction problem via gauge-fixing and rotational equivariance. The mathematical content is substantial and the analytical pipeline is sound. However, the document is currently in a **research notes / working report** stage rather than a **formal research article**. The gap between the two is significant but bridgeable. The main issues fall into six categories: (1) incomplete mathematical rigor, (2) structural/organizational problems, (3) missing standard article components, (4) figure and presentation issues, (5) LaTeX/typesetting problems, and (6) mathematical content gaps and potential errors.

---

## 1. Critical Mathematical Issues

### 1.1 Incomplete Proofs and Unverified Claims

**Priority: HIGHEST**

The paper's central argument chain has three explicitly stated open problems (Section 5, lines 1004–1010):

1. **Rotational equivariance of $F_i$** — The claim that $F_i(R_\omega \delta y, R) = R_\omega F_i(\delta y, R)$ (line 949) is asserted but never proved. This is the linchpin of the entire kernel resolution strategy. Without it, the gauge-fixing argument collapses. A formal proof should verify this identity by:
   - Writing $F_i$ explicitly in terms of the rotation-transformed variables.
   - Showing the Hessian $H_V$ (diagonal in the eigenvector basis) commutes appropriately with rotations when composed with the $Q_R$ interaction terms.
   - Handling the $R > 0$ terms carefully, since $Q_R$ depends on $|\tilde{y}_k - \tilde{y}_i|$ which is rotation-invariant, but the directional dependence $(y_k - y_i)/|y_k - y_i|$ must also transform correctly.

2. **Invertibility of the linearized operator on the gauge-fixed subspace** — This requires showing that after projecting out the translational and rotational kernel directions, $D_{\delta y}F_i(\delta y^0, 0)$ restricted to $\{w_{\text{rot}}\}^\perp \cap \Omega^\perp$ is nonsingular. This is a finite-dimensional linear algebra problem (a $4 \times 4$ matrix after gauge-fixing from $6 \times 6$) and should be computable explicitly for the equilateral case.

3. **Extension to continuous $\theta$** — The current treatment handles the aligned equilateral case ($\theta = 0$) in the perturbation analysis. The $\theta$-dependent case requires showing the IFT applies uniformly over the admissible $\theta$-range.

**Recommendation:** These three items must be resolved or clearly delineated as the paper's stated open problems with a precise "Conjecture + partial evidence" framing, or completed in full. A publication in a strong journal will require at least items 1 and 2 to be complete.

### 1.2 Claim/Theorem Hierarchy Issues

**Priority: HIGH**

The document uses `Claim` environments for statements of varying importance and rigor:

- **Claim 1** (line 241): "$x_0$ is a non-degenerate critical point of $V$" — This is an *assumption*, not a claim. It should be a **Hypothesis** or **Assumption** environment.
- **Claim 2** (line 246): "$Rz_i(R) \to x_0$" — This is a known result from Kirr–Natarajan [2018]. It should be a **Proposition** or **Theorem** with proper attribution.
- **Claim 3** (line 722): The explicit solution $\delta y^0$ — This is a genuine result with a proof provided. It should be a **Proposition** or **Lemma**.
- **Claim 4** (line 842): Nontrivial kernel — Another genuine result with proof. Should be a **Proposition**.
- **Claim 5** (line 920): Gauge-fixing yields solutions — This is stated as a claim but the "proof" (lines 926–945) actually only shows that $\sum_i F_i = 0$ recovers the reduced system, not that the gauge-fixing constraints yield unique solutions. This needs significant strengthening.

**Recommendation:** Adopt the standard hierarchy: **Assumption** → **Definition** → **Lemma** → **Proposition** → **Theorem** → **Corollary**. Reserve `Claim` for statements that are conjectured but not yet proved.

### 1.3 Potential Mathematical Error: Sign Convention for $\chi$

**Priority: HIGH**

In the interaction coefficient definitions, the paper defines (line 740):
$$\chi = \max\{\alpha_{12}, \alpha_{13}, \alpha_{23}\} = \alpha_{12} = \alpha_{13}$$

But earlier (line 1127), $\chi$ is defined as:
$$\chi = \frac{2Q_R(m_R)}{R^4 m_R}$$

The relationship between these two definitions needs to be explicitly established. The normalization convention $\alpha_{ki} = Q_R(|z_k - z_i|)/Q_R(m_R)$ means $\max \alpha_{ki} = 1$ (attained at the minimum distance), so $\chi$ as defined in the notation table should equal $\alpha_{12} \cdot \chi_{\text{table}}$, not $\alpha_{12}$ itself. This conflation is confusing and potentially erroneous. The proof of Claim 3 appears to use $\chi = \alpha_{12}$ (the interaction coefficient value) rather than the overall interaction strength. Clarify.

### 1.4 The $\tilde{\alpha}_{ki}$ vs $\alpha_{ki}$ Confusion

**Priority: HIGH**

The paper defines two related quantities:
- $\alpha_{ki} = Q_R(|z_k - z_i|)/Q_R(m_R)$ — the ratio of interaction functionals.
- $\tilde{\alpha}_{ki} = e^{-\langle y_k - y_i, \delta y_k - \delta y_i \rangle}$ — the perturbed interaction coefficient at $R = 0$.

In the proof of Claim 3 (line 737–740), it is shown that $\tilde{\alpha}_{12} = \tilde{\alpha}_{13} = 1$ and $\tilde{\alpha}_{23} = \alpha_{23}/\alpha_{12}$. The notation switches between these without clear flagging. In the $F_i(\delta y, 0)$ system (eq. 11), the coefficient is $\tilde{\alpha}_{ki}\chi$, but in the algebraic system being verified, the coefficients used are $\alpha_{12}$ and $\alpha_{23}/\alpha_{12}$. The relationship $\tilde{\alpha}_{ki} \cdot \chi = \alpha_{ki}$ (where $\chi = \alpha_{12}$) must be stated explicitly.

### 1.5 Missing Decay Conditions on $V$

**Priority: MEDIUM-HIGH**

Line 235–238 contains:
```latex
Under certain decay conditions on the potential $V$ as $|x| \to \infty$
\begin{comment}
List them
\end{comment}
```

This is a commented-out TODO that made it into the document. The decay conditions on $V$ are essential hypotheses — without them, the concentration compactness argument and the boundedness of $Rz_i$ are not justified. These conditions should include at minimum:
- $V \in L^\infty(\mathbb{R}^n)$ or $V \in L^p + L^\infty$ for appropriate $p$.
- $V(x) \to 0$ as $|x| \to \infty$ (or at least $V$ bounded and approaching a limit).
- The non-degeneracy hypothesis on $x_0$.

### 1.6 The Corollary About the Unstable Manifold

**Priority: MEDIUM**

The Corollary on line 326–329 states that $\text{span}\{y_1, \ldots, y_M\} \subseteq U$ where $U$ is the span of eigenvectors with negative eigenvalues. The notation says "unstable manifold" but this is actually the **unstable subspace** of the linearized system at $x_0$ (which is a linear subspace, not a manifold). Furthermore, the result follows from the algebraic system $H_V y_i = \sum \alpha_{ki}\chi(y_k - y_i)$ by projecting onto eigendirections of $H_V$, but this argument is not written out.

---

## 2. Structural and Organizational Issues

### 2.1 Section Structure

**Priority: HIGH**

The current structure is:

1. Introduction
2. Background
3. Summary of Scenarios
4. Research Problem
5. Conclusion
6. Notation

**Problems:**
- "Summary of Scenarios" is not a standard section name. This is really "Classification of Multi-Peak Configurations" or "Admissible Geometries."
- "Research Problem" conflates the formulation of the perturbation system with its analysis and partial resolution.
- The Notation section (currently Section 6) containing 10 tables spanning ~250 lines is excessive for a journal article. Some journals may accept a notation table in an appendix, but most expect notation to be introduced inline.
- There is no "Main Results" theorem statement.

**Recommended structure for a journal article:**

1. **Introduction** (motivation, informal statement of results, literature context)
2. **Setup and Preliminaries** (NLS equation, rescaling, concentration compactness decomposition, Lyapunov–Schmidt reduction — all attributed to prior work)
3. **Main Results** (formal theorem/proposition statements)
4. **Classification of Limiting Configurations** ($R = 0$ analysis: collinear, isosceles, equilateral, rotational)
5. **The Perturbation System** (formulation of $F_i$, explicit solution at $R = 0$)
6. **Kernel Analysis and Resolution** (kernel structure, gauge-fixing, projections, equivariance)
7. **Conclusion and Open Problems**
- **Appendix A:** Notation Reference (if desired)
- **Appendix B:** Detailed Computations (the long matrix verifications)

### 2.2 Missing Literature Review / Related Work

**Priority: HIGH**

The introduction mentions Lions (1984), Kirr–Natarajan (2018), and briefly Floer–Weinstein (1986), but does not discuss:

- **Multi-bump solutions literature:** The work of Séré (1992), Coti Zelati–Rabinowitz (1992), and subsequent multi-bump/multi-peak constructions for NLS and related equations. The Kirr–Knox work should be positioned relative to these.
- **Lyapunov–Schmidt reduction for NLS:** Wei (1996, 1997), Ambrosetti–Malchiodi–Ni (2004) and the extensive literature on spike solutions of singularly perturbed NLS. How does the Kirr–Natarajan framework compare?
- **Concentration on potential wells:** Byeon–Wang (2002, 2005), del Pino–Felmer (1996, 1997) — especially their results on the number and location of peaks.
- **The optical fiber context:** Akhmediev–Ankiewicz (1997), Agrawal's textbook on nonlinear fiber optics. The physical motivation section should cite these standard references.

**Recommendation:** Add a paragraph in the Introduction comparing and contrasting with the multi-bump literature. Emphasize what is new: the explicit classification of admissible peak *geometries* (not just existence), the interaction between Hessian eigenvalue ratios and peak arrangement, and the systematic kernel decomposition.

### 2.3 Abstract Needs Tightening

**Priority: MEDIUM**

The abstract (lines 61–91) is 267 words, which is at the upper end for most journals. More importantly, it:
- Spends significant space on methodology description rather than stating results.
- Does not clearly state the main theorem(s) in a self-contained way.
- Uses "$\tfrac{1}{3} \leq \mu \leq 3$" which is a specific result that should be contextualized.

**Recommendation:** Restructure as: (1) Problem statement (2 sentences), (2) Main results (4–5 sentences with the key classification theorem stated), (3) Method/approach (2 sentences). Target 200 words.

---

## 3. Missing Standard Article Components

### 3.1 No Formal Theorem Statements

**Priority: HIGHEST**

A mathematics research article requires clearly stated, numbered theorems with precise hypotheses. The current report has no `\newtheorem{thm}{Theorem}` environment and no theorem statements. The main results should be consolidated into:

**Theorem 1** (Classification of limiting configurations): *Under Hypotheses (H1)–(H3), the admissible $M$-peak configurations at a non-degenerate critical point $x_0$ of $V$ with Hessian eigenvalues $\lambda_1, \lambda_2 < 0$ are:*
- *(a) Collinear: ...*
- *(b) Isosceles triangle ($M=3$): requires $\lambda_2 = 3\lambda_1$...*
- *(c) Equilateral triangle ($M=3$): requires $|\lambda_2|/|\lambda_1| \geq 3$...*
- *(d) Rotational equilateral ($M=3$): admissible $(\mu, \theta)$ regions given by Table 1...*

**Theorem 2** (Perturbation to finite energy): *For the equilateral triangle configuration with $\mu > 3$, there exists $\varepsilon > 0$ and a smooth family $\delta y(R)$ for $R \in [0, \varepsilon)$ solving $F_i(\delta y, R) = 0$ with $\delta y(0) = \delta y^0$, provided the rotational equivariance condition holds.*

### 3.2 No Keywords or MSC Codes

**Priority: MEDIUM**

Every journal submission needs:
- **Keywords:** nonlinear Schrödinger equation, multi-peak solutions, concentration compactness, Lyapunov–Schmidt reduction, ground states, optical fibers
- **MSC 2020:** 35Q55 (NLS equations), 35B40 (Asymptotic behavior), 35B32 (Bifurcations), 37K50 (Bifurcation problems), 78A60 (Lasers, masers, optical bistability)

### 3.3 No Acknowledgments

**Priority: LOW**

Standard practice: acknowledge funding sources, helpful discussions, etc.

### 3.4 No Running Headers/Footers

**Priority: LOW**

The `extarticle` class is fine for a report but most journals provide their own class files. For arXiv preprint submission, `amsart` or a standard `article` class with `amsmath` is more conventional.

---

## 4. Figure and Presentation Issues

### 4.1 Hardcoded Local File Paths

**Priority: HIGHEST (blocks compilation)**

Every `\includegraphics` command uses an absolute path to a local SMB share:
```
/run/user/1000/gvfs/smb-share:server=skynet.local,share=public/nls-research/figures/output/
```

This means the document **will not compile on any other machine**. Replace all paths with relative paths:
```latex
\includegraphics[width=\linewidth]{figures/figure_7_methodology_flowchart}
```

### 4.2 Missing/Placeholder Figure Captions

**Priority: HIGH**

- Figure 2.2 (line 253): Caption is literally `"graphic"`. This needs a real caption describing the concentration compactness visualization.
- Several other captions are too terse: "Peak Configurations" (Figure 3.1), "$\alpha$ vs. $\theta$ plots for various $\mu$" (Figure 3.4). Captions should be self-contained mini-descriptions.

### 4.3 Figure Quality and Consistency

**Priority: MEDIUM**

- All figures use `[width=1\linewidth]` or `[width=0.75\linewidth]` with `[H]` (forced here) placement. For a journal article, use `[htbp]` and let LaTeX handle placement.
- The pipeline flowchart (Figure 1.1) is placed at the end of the Introduction but should arguably appear earlier or in a methods/overview section.
- Figures are PNG/PDF from matplotlib. For journal submission, EPS or high-resolution PDF vector graphics are preferred. The generation scripts already produce these formats, so just ensure the LaTeX points to the vector versions.

### 4.4 No Figure Cross-References

**Priority: MEDIUM**

The text does not reference most figures with `\ref{}`. Figures should be referenced in the text near where they appear: "As shown in Figure~\ref{fig:admissibility}, the admissible region..."

---

## 5. LaTeX and Typesetting Issues

### 5.1 Character Encoding

**Priority: HIGH**

The document uses `\xf6` (raw byte for ö) instead of `\"o` or proper UTF-8 encoding. This is a LyX export artifact. All instances of `\xf6` should be replaced with `\"o` (or switch to `\usepackage[utf8]{inputenc}`). Affected instances: "Schrödinger" appears ~10 times.

Similarly `\xe9` for é in "Fréchet" (line 1269).

### 5.2 Equation Numbering

**Priority: MEDIUM**

The first equation is manually tagged `\tag{1}\label{eq:NLS}` instead of using automatic numbering. Later equations use `\label` without `\tag`, creating inconsistency. Remove all `\tag` commands and let LaTeX number equations automatically via `equation` environments.

### 5.3 Use of `\phantom{}`

**Priority: LOW-MEDIUM**

Lines 245, 284 use `\phantom{}` to create vertical spacing between theorem environments. This is a hack. Use `\medskip`, or better, define proper spacing in the theorem environments.

### 5.4 The `\Downarrow` Transitions

**Priority: MEDIUM**

The document uses standalone display-math `\Downarrow` arrows to indicate logical flow (lines 390, 404, 460, 527, etc.). While common in lecture notes, this is not standard in published mathematics. Replace with prose transitions: "which yields," "it follows that," "combining these we obtain," etc.

### 5.5 Missing `\qedhere` and Proof Formatting

**Priority: LOW**

Some proofs end abruptly. Ensure all `proof` environments end with a clear QED symbol, and that the final line of each proof reads naturally.

### 5.6 Document Class

**Priority: MEDIUM**

`extarticle` is nonstandard for math publications. Switch to `amsart` (AMS article class) which provides:
- Standard theorem environments (`\newtheorem`)
- Proper AMS math formatting
- Standard section numbering
- Built-in bibliography support

---

## 6. Mathematical Content to Add or Strengthen

### 6.1 The Collinear Case Needs More Development

**Priority: MEDIUM**

The collinear case (Example 3.1, lines 340–356) is treated very briefly with a generic formula. For a complete classification paper, this should include:
- Explicit construction for $M = 2$ (two-peak, which is the simplest and most important case).
- The chain of interaction coefficients $\alpha_{i,i+1}$ and the resulting constraints.
- Connection to the known results for symmetric double-well potentials.

### 6.2 The Isosceles Case Is Incomplete

**Priority: MEDIUM**

Example 3.2 derives $\alpha_{12} = \alpha_{13} = -\lambda_2/3$ and $\lambda_2 = 3\lambda_1$ but does not:
- Discuss what happens when this eigenvalue ratio is not met (the configuration is inadmissible).
- Explain the geometric meaning: the isosceles triangle degenerates to equilateral when $2a = 1$, so the isosceles case with $2a > 1$ is a strictly non-equilateral triangle.
- Verify the binding constraint $\alpha_{23}$ (since $y_2$ and $y_3$ are not nearest neighbors in the isosceles case).

### 6.3 Stability Discussion Is Too Brief

**Priority: MEDIUM**

The Conclusion (line 1013–1018) mentions that multi-peak branches have "at least two negative directions in the linearized operator, implying orbital instability." This is a significant statement that needs:
- A precise reference to which linearized operator (the full PDE linearization at $\psi_E$).
- Citation to GSS (1987, 1990) for the abstract stability/instability criterion.
- Discussion of whether the instability is generic or depends on the configuration type.

### 6.4 Missing: Physical Interpretation Section

**Priority: MEDIUM**

The Introduction motivates the problem via optical fibers, but the paper never returns to physics. A brief discussion section should interpret the results physically:
- What does $\mu = |\lambda_2|/|\lambda_1|$ mean for fiber design?
- What refractive-index profiles correspond to the different configuration types?
- Are the multi-peak states experimentally observable, given their instability?

---

## 7. Repository and Reproducibility

### 7.1 Build System

**Priority: MEDIUM**

Add a `Makefile` or `latexmk` configuration that:
1. Runs the Python figure scripts.
2. Compiles the LaTeX document with proper bibtex passes.
3. Is documented in the README.

### 7.2 Mathematica Notebooks

**Priority: MEDIUM**

The existing notebooks (`Research-Triangles.nb`, `Research-Triangle-Pertubation.nb`) are minimal — the perturbation notebook contains a single cell that verifies `True` via `Reduce`. These should be expanded into a comprehensive verification suite (see the audit notebook delivered with this report).

### 7.3 Version Control

**Priority: LOW**

The repo contains LyX backup files (`Research.lyx~`, `#Research.lyx#`) that should be in `.gitignore`. The `archive/` directory with old drafts is fine for internal use but should be excluded from any public/submission version.

---

## 8. Prioritized Action Plan

### Phase 1: Mathematical Completion (Weeks 1–4)
1. Prove rotational equivariance of $F_i$ or state it precisely as a conjecture.
2. Compute the linearized operator on the gauge-fixed subspace explicitly and verify invertibility.
3. Resolve all Claim → Theorem/Proposition/Assumption conversions.
4. Fill in the missing decay conditions on $V$.
5. Clarify the $\chi$ / $\alpha_{ki}$ / $\tilde{\alpha}_{ki}$ normalization.

### Phase 2: Structural Rewrite (Weeks 3–5)
1. Reorganize into the recommended section structure.
2. Write formal theorem statements with precise hypotheses.
3. Add literature review paragraph.
4. Add keywords, MSC codes.
5. Tighten the abstract.

### Phase 3: Presentation Polish (Weeks 5–6)
1. Fix all figure paths to relative.
2. Write proper captions for all figures.
3. Replace `\Downarrow` transitions with prose.
4. Fix encoding issues (`\xf6` → `\"o`).
5. Switch to `amsart` or target journal class.
6. Move notation tables to appendix; introduce notation inline.
7. Add physical interpretation discussion.

### Phase 4: Final Preparation (Week 7)
1. Run the Mathematica audit notebook and verify all computations.
2. Proofread for typos and mathematical consistency.
3. Ensure all figures are vector format.
4. Prepare arXiv submission package (single `.tex` file with all figures).
5. Write cover letter for journal submission.

---

## 9. Target Venues (Ranked by Fit)

1. **Journal of Differential Equations** — Good fit for classification + existence results for NLS.
2. **Nonlinearity** (IOP) — Strong fit for the geometric classification angle.
3. **Communications in Mathematical Physics** — If the full existence theorem (including kernel resolution) is completed.
4. **SIAM Journal on Mathematical Analysis** — If the optical fiber application is emphasized.
5. **arXiv preprint** (math.AP) — Recommended as a first step regardless of journal choice.

---

## 10. Summary of Critical Items

| # | Issue | Priority | Effort |
|---|-------|----------|--------|
| 1 | Prove rotational equivariance | HIGHEST | High |
| 2 | Formal theorem statements | HIGHEST | Medium |
| 3 | Fix figure paths (blocks compilation) | HIGHEST | Low |
| 4 | Prove invertibility on gauge-fixed subspace | HIGH | High |
| 5 | Resolve Claim hierarchy | HIGH | Low |
| 6 | Clarify $\chi$ normalization | HIGH | Low |
| 7 | Fill in $V$ decay conditions | HIGH | Low |
| 8 | Add literature context | HIGH | Medium |
| 9 | Restructure sections | HIGH | Medium |
| 10 | Fix `\xf6` encoding | HIGH | Low |
| 11 | Write proper figure captions | HIGH | Low |
| 12 | Expand collinear/isosceles cases | MEDIUM | Medium |
| 13 | Add physical interpretation | MEDIUM | Medium |
| 14 | Switch to `amsart` class | MEDIUM | Low |
| 15 | Move notation to appendix | MEDIUM | Low |
