# Multi-Peak Ground States of the Nonlinear Schrödinger Equation in Optical Fibers

Mathematical research on multi-peak ground-state solutions of the stationary NLS with external potential in the high-energy regime ($E \to \infty$), using concentration compactness and Lyapunov-Schmidt reduction. The work establishes existence, uniqueness, and stability of ground-state configurations with multiple concentration peaks at a nondegenerate minimum of the potential $V$.

## Authors

- **Eduard Kirr, PhD** — Associate Professor, University of Illinois at Urbana-Champaign
- **Ethan Knox, BS** — Master's Student, University of Illinois at Urbana-Champaign

## Overview

We study the stationary NLS equation with external potential in the high-energy regime:

$$(-\Delta + V + E)\psi + \sigma|\psi|^{2p}\psi = 0$$

as $E \to \infty$. The work is motivated by optical fiber physics — the envelope equation in silica fibers with Kerr nonlinearity maps to this stationary Schrödinger equation, where the refractive-index profile becomes an external potential $V$.

### Key Framework

1. **Rescaling.** Define $u_E(x) = E^{-1/(2p)} \psi_E(E^{-1/2} x)$, which transforms the problem into a family of equations parametrized by $R = E^{-1/2} \to 0$.

2. **Concentration compactness decomposition.** Each solution decomposes into $M$ well-separated profiles $Q$ centered at distinct points $z_1, \ldots, z_M$, with a remainder $h_R$:

$$\psi_E(x) = \sum_{i=1}^{M} u_E(x - Rz_i) + h_R(x)$$

3. **Lyapunov-Schmidt reduction.** Projection onto the kernel of the linearized operator yields a finite-dimensional algebraic system for peak positions:

$$H_V(x_0) y_i = \sum_{k \in N_i} \alpha_{ki} \chi(y_k - y_i)$$

where $H_V(x_0)$ is the Hessian of $V$ at the critical point $x_0$, and $\alpha_{ki}$ are interaction coefficients between peaks $k$ and $i$.

4. **Perturbation system.** The full problem at finite energy $R > 0$ is captured by $F_i(\delta y, R) = 0$, with $\delta y$ denoting peak displacements from the $R = 0$ limiting configuration.

## Main Results

### 1. Classification of Limiting Configurations ($R = 0$)

The algebraic system admits four qualitatively distinct configuration types:

- **Collinear** — Peaks align along a single eigenvector of the Hessian $H_V(x_0)$. Applies for $M = 2$ and $M \geq 3$.
- **Isosceles triangle** ($M = 3$) — Two equal interaction lengths; requires $\lambda_2 = 3\lambda_1$ (eigenvalue ratio constraint).
- **Equilateral triangle** ($M = 3$) — All three peaks equidistant; requires $|\lambda_2|/|\lambda_1| \geq 3$.
- **Equilateral with rotational degree of freedom** ($M = 3$) — Admissibility regions in the $(\mu, \theta)$ parameter space, where $\mu = |\lambda_2|/|\lambda_1|$ and $\theta$ is the orientation angle.

**Admissibility summary:**

| $\mu$-range | $\theta$-range | Configuration |
|---|---|---|
| $\mu = 1$ | All $\theta \in [0, 2\pi)$ | Fully admissible (isotropic potential) |
| $1/3 \leq \mu \leq 3$ | All $\theta \in [0, 2\pi)$ | All orientations admissible |
| $\mu < 1/3$ or $\mu > 3$ | Restricted angular windows | Partially admissible |

The admissible $(\mu, \theta)$ phase diagram is shown in Figure 1.

### 2. Perturbation System ($R > 0$)

- **Formulation.** The perturbation system $F_i(\delta y, R) = 0$ extends the $R = 0$ algebraic system to finite energy, with explicit continuation from $R = 0$ via implicit function theorem arguments.
- **Explicit solution at $R = 0$.** For the equilateral triangle configuration, the leading-order displacement corrections are:

$$\delta y_1^0 = -\gamma r_2, \qquad \delta y_{2,3}^0 = \mp\beta r_1 + \frac{\gamma}{2}r_2$$

where $\beta = \frac{1}{2}\ln(\alpha_{12}/\alpha_{23})$ and $\gamma = 2\beta/(3\sqrt{3})$.

- **Nontrivial kernel.** The linearized operator has a kernel $\Omega \supseteq \{(v, v, v) : v \in \mathbb{R}^2\}$ arising from translational invariance of the equation.

### 3. Kernel Resolution Strategy

- **Projection decomposition.** The perturbation system is decomposed into $P_\perp F$ (transverse projection) and $P_\parallel F$ (kernel projection) components.
- **Gauge-fixing constraints.** The kernel is resolved by imposing $\delta y_1 \perp r_1$ and $\delta y_3 - \delta y_2 \parallel r_1$, which remove the translational zero modes.
- **Rotational equivariance.** The structural property enabling application of the implicit function theorem on the gauge-fixed complement.

### Open Tasks

1. Verify rotational equivariance of $F_i$ rigorously.
2. Prove invertibility of the linearized operator on the gauge-fixed subspace.
3. Extend to configurations with continuous rotation angle $\theta$.

## Repository Structure

```
nls-research/
├── report/                              # Main research paper
│   ├── Research.lyx                     # Primary LyX source (5446 lines)
│   ├── Research.tex                     # Exported LaTeX source
│   ├── Research.pdf                     # Compiled PDF
│   ├── Research.bib                     # BibTeX bibliography (14 entries)
│   ├── Research-old.tex                 # Previous version (LaTeX)
│   ├── Research-old.pdf                 # Previous version (PDF)
│   └── archive/                         # Archived drafts
│       ├── Research_Draft_2.3.lyx
│       ├── Research_Draft_2.3.tex
│       ├── Research_Draft_2.3.pdf
│       ├── mu-theta-feasibility-plot-notitle.eps
│       ├── mu-theta-feasibility-plot-notitle.png
│       └── mu-theta-feasibility-plot-title.png
│
├── figures/                             # Figure generation (Python)
│   ├── config.py                        # Output directory configuration
│   ├── main.py                          # Entry point: generates all 7 figures
│   ├── output/                          # Generated outputs (EPS, PDF, PNG)
│   └── figures/                         # Individual figure modules
│       ├── figure_1.py                  # Admissibility & binding constraints (μ-θ phase diagram)
│       ├── figure_2.py                  # Peak configuration diagrams
│       ├── figure_3.py                  # α vs θ plots for various μ
│       ├── figure_4.py                  # Multi-peak profile decomposition
│       ├── figure_5.py                  # Concentration compactness visualization
│       ├── figure_6.py                  # Interaction functional Q_R decay
│       ├── figure_7.py                  # Analytical pipeline flowchart
│       └── figure_7.md                  # Markdown notes for figure 7
│
├── triangle/                            # Mathematica computations
│   ├── Research-Triangles.nb            # Triangle configuration analysis
│   ├── Research-Triangle-Pertubation.nb # Perturbation analysis
│   ├── Research-Triangle-Pertubation.pdf # Compiled perturbation notebook
│   └── triangle_soln.md                 # Markdown derivation of admissibility
│
├── eddy-notes/                          # Research notes (gitignored)
│   ├── ResearchNotes.md
│   ├── Peak interaction project.pdf
│   └── Multipeaks at same critical point.pdf
│
├── reference/                           # Reference papers (14 PDFs)
│
├── requirements.txt                     # Python dependencies
├── .forgejo/workflows/                  # CI: GitHub mirror
└── .gitignore
```

## Figures

| # | Filename | Description |
|---|----------|-------------|
| 1 | `figure_1_admissibility_binding` | $(\mu, \theta)$ phase diagram showing admissible regions and binding constraints |
| 2 | `figure_2_peak_configurations` | Peak arrangements: collinear ($M=2$, $M=5$), isosceles, equilateral, perturbation |
| 3 | `figure_3_alpha_vs_theta_vs_mu` | Interaction coefficients $\alpha_{12}$, $\alpha_{13}$, $\alpha_{23}$ vs $\theta$ for $\mu \in \{0.3, 1.0, 3.0, 6.0\}$ |
| 4 | `figure_4_profile_decomposition` | Multi-peak profile decomposition with remainder $h_R$ |
| 5 | `figure_5_concentration_vs_peaks` | Concentration compactness: $Rz_i \to x_0$ as $E \to \infty$ |
| 6 | `figure_6_interaction_decay` | Interaction functional $Q_R$ exponential decay and $\alpha_{ki}$ coefficient |
| 7 | `figure_7_methodology_flowchart` | Analytical pipeline: Kirr-Natarajan (2018) → Kirr-Knox (2026) |

All figures are generated in EPS, PDF, and PNG formats in `figures/output/`.

## Prerequisites

### Python

- Python 3.8+
- numpy
- matplotlib

### LaTeX / LyX

- LyX 2.4+ (for editing `.lyx` source)
- A standard LaTeX distribution (TeX Live or MiKTeX) with packages: `extarticle`, `geometry`, `amsmath`, `amsthm`, `amssymb`, `graphicx`, `booktabs`, `mathtools`, `varwidth`, `units`, `float`, `verbatim`

### Mathematica

- Wolfram Mathematica (for `triangle/*.nb` notebooks)

## Setup and Usage

### Generating Figures

```bash
cd figures
pip install -r ../requirements.txt
python main.py
```

Figures are written to `figures/output/` in EPS, PDF, and PNG formats.

### Compiling the Paper

From `report/`:

```bash
pdflatex Research.tex
bibtex Research
pdflatex Research.tex
pdflatex Research.tex
```

Or open `Research.lyx` in LyX and use the built-in export/compile.

## Background and Context

### Physical Motivation

The propagation of intense light through optical fibers is governed by the interplay between linear dispersion and Kerr nonlinearity ($\chi^{(3)}$). Starting from Maxwell's equations, the slowly varying envelope approximation yields the NLS equation:

$$i\partial_z A - \frac{\beta_2}{2}\partial_T^2 A + \gamma |A|^2 A = 0$$

Spatial variations in the refractive index (fiber Bragg gratings, graded-index fibers, photonic-crystal fibers) enter as an external potential $V(\mathbf{x})$, making the transverse modal equation mathematically equivalent to the stationary Schrödinger equation.

### Mathematical Framework

The analysis builds on:

- **Concentration compactness** (Lions, 1984) — classifying minimizing sequences in variational problems on unbounded domains.
- **Global bifurcation analysis** (Kirr and Natarajan, 2018) — rescaling and multi-peak decomposition of ground states via equivariant bifurcation theory.
- **Lyapunov-Schmidt reduction** — reduction to a finite-dimensional algebraic system governing peak positions.
- **Grillakis-Shatah-Strauss stability theory** — orbital instability of multi-peak states in Hamiltonian systems with symmetry.

## References

1. [Hasegawa and Tappert, 1973a] A. Hasegawa and F. Tappert. Transmission of stationary nonlinear optical pulses in dispersive dielectric fibers. I. Anomalous dispersion. *Applied Physics Letters*, 23(3):142–144, 1973.

2. [Hasegawa and Tappert, 1973b] A. Hasegawa and F. Tappert. Transmission of stationary nonlinear optical pulses in dispersive dielectric fibers. II. Normal dispersion. *Applied Physics Letters*, 23(4):171–172, 1973.

3. [Lions, 1984a] P. L. Lions. The concentration-compactness principle in the Calculus of Variations. The locally compact case, part 1. *Annales de l'Institut Henri Poincaré. Analyse non linéaire*, 1(2):109–145, 1984.

4. [Lions, 1984b] P. L. Lions. The concentration-compactness principle in the Calculus of Variations. The locally compact case, part 2. *Annales de l'Institut Henri Poincaré. Analyse non linéaire*, 1(4):223–283, 1984.

5. [Kirr and Natarajan, 2018] E. Kirr and V. Natarajan. The global bifurcation picture for ground states in nonlinear Schrödinger equations. arXiv preprint arXiv:1811.05716, 2018.

6. [Berestycki and Lions, 1983] H. Berestycki and P. L. Lions. Nonlinear scalar field equations, I existence of a ground state. *Archive for Rational Mechanics and Analysis*, 82(4):313–345, 1983.

7. [Floer and Weinstein, 1986] A. Floer and A. Weinstein. Nonspreading wave packets for the cubic Schrödinger equation with a bounded potential. *Journal of Functional Analysis*, 69(3):397–408, 1986.

8. [Kirr, Kevrekidis, and Pelinovsky, 2011] E. Kirr, P. G. Kevrekidis, and D. E. Pelinovsky. Symmetry-breaking bifurcation in the nonlinear Schrödinger equation with symmetric potentials. *Communications in Mathematical Physics*, 308(3):795–844, 2011.

9. [Kirr, 2016] E. Kirr. Long time dynamics and coherent states in nonlinear wave equations. arXiv preprint arXiv:1605.08167, 2016.

10. [Grillakis, Shatah, and Strauss, 1987] M. Grillakis, J. Shatah, and W. Strauss. Stability theory of solitary waves in the presence of symmetry, I. *Journal of Functional Analysis*, 74(1):160–197, 1987.

11. [Grillakis, Shatah, and Strauss, 1990] M. Grillakis, J. Shatah, and W. Strauss. Stability theory of solitary waves in the presence of symmetry, II. *Journal of Functional Analysis*, 94(2):308–348, 1990.

12. [Kirr and Zarnescu, 2007] E. Kirr and A. Zarnescu. On the asymptotic stability of bound states in 2D cubic Schrödinger equation. *Communications in Mathematical Physics*, 272(2):443–468, 2007.

13. [Kirr and Zarnescu, 2009] E. Kirr and A. Zarnescu. Asymptotic stability of ground states in 2D nonlinear Schrödinger equation including subcritical cases. arXiv preprint arXiv:0805.3888, 2009.

14. [Kirr and Mizrak, 2008] E. Kirr and O. Mizrak. Asymptotic stability of ground states in 3D nonlinear Schrödinger equation including subcritical cases. arXiv preprint arXiv:0803.3377, 2008.

## License

This project is for academic research purposes. Contact the authors for licensing information.
