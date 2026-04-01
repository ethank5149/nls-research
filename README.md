# Multi-Peak Ground States of the Nonlinear Schrödinger Equation with External Potential

Mathematical research on multi-peak ground-state solutions of the stationary NLS with external potential in the high-energy regime ($E \to \infty$), using concentration compactness and Lyapunov-Schmidt reduction. The work establishes existence and classification of ground-state configurations with multiple concentration peaks at a non-degenerate local maximum of the potential $V$.

## Authors

- **Eduard Kirr** — Associate Professor, University of Illinois at Urbana-Champaign
- **Ethan Knox** — Graduate Student, University of Illinois at Urbana-Champaign

## Overview

We study the stationary NLS equation with external potential in the high-energy regime:

$$(-\Delta + V + E)\psi + \sigma|\psi|^{2p}\psi = 0$$

as $E \to \infty$. The work is motivated by optical fiber physics — the envelope equation in silica fibers with Kerr nonlinearity ($\chi^{(3)}$) maps to this stationary Schrödinger equation, where the refractive-index profile becomes an external potential $V(\mathbf{x})$.

### Hypotheses on the Potential

Following Kirr, Kevrekidis, and Pelinovsky (2011), Floer and Weinstein (1986), and Kirr and Natarajan (2018):

- **(V1)** $V \in L^\infty(\mathbb{R}^n)$ — bounded and measurable.
- **(V2)** $\lim_{|x|\to\infty} V(x) = 0$ — vanishes at infinity.
- **(V3)** $V$ is $C^2$ near $x_0$, a non-degenerate critical point: $\nabla V(x_0) = 0$ and $H_V(x_0)$ invertible.
- **(V4)** All eigenvalues of $H_V(x_0)$ are strictly negative — $x_0$ is a strict local maximum.

### Key Framework

1. **Rescaling.** Define $u_E(x) = E^{-1/(2p)} \psi_E(E^{-1/2} x)$, which transforms the problem into a family of equations parametrized by $R = E^{-1/2} \to 0$.

2. **Concentration compactness decomposition.** Each solution decomposes into $M$ well-separated profiles $u_R$ centered at distinct points $z_1(R), \ldots, z_M(R)$, with a remainder $h_R$:

$$u_E(x) = \sum_{i=1}^{M} u_R(x - z_i(R)) + h_R(x)$$

3. **Lyapunov-Schmidt reduction.** Projection onto the kernel of the linearized operator yields a finite-dimensional algebraic system for peak positions:

$$H_V(x_0)\, y_i = \chi \sum_{k \in N_i} \alpha_{ki}(y_k - y_i)$$

where $H_V(x_0)$ is the Hessian of $V$ at the critical point $x_0$, $\alpha_{ki} = Q_R(|z_k - z_i|)/Q_R(m_R) \in [0,1]$ are normalized interaction coefficients, and $\chi > 0$ is the overall interaction strength.

4. **Perturbation system.** The full problem at finite energy $R > 0$ is captured by $F_i(\delta y, R) = 0$, with $\delta y$ denoting peak displacements from the $R = 0$ limiting configuration.

## Main Results

### 1. Classification of Limiting Configurations ($R = 0$)

**Theorem 3.1.** Under Hypotheses (V1)–(V4), the admissible solutions of the algebraic system with all $\alpha_{ki} \geq 0$ are:

- **(a) Collinear** ($M \geq 2$) — Peaks aligned along a single eigenvector $\mathbf{r}_1$ of $H_V(x_0)$.
- **(b) Isosceles triangle** ($M = 3$) — Peaks with $\alpha_{12} = \alpha_{13} = -\lambda_2/3$; requires the eigenvalue constraint $\lambda_2 = 3\lambda_1$.
- **(c) Equilateral triangle** ($M = 3$) — All peaks at unit mutual distance; $\alpha_{12} = \alpha_{13} = -\lambda_2/3$, $\alpha_{23} = \lambda_2/6 - \lambda_1/2$; requires $|\lambda_2|/|\lambda_1| \geq 3$.
- **(d) Rotational equilateral** ($M = 3$, angle $\theta$) — Equilateral triangle rotated by $\theta$ relative to the eigenvector basis, with admissible $(\mu, \theta)$ regions given below.

**Admissibility regions** (where $\mu = |\lambda_2|/|\lambda_1|$):

| $\mu$-Range | $\theta$-Range |
|---|---|
| $0 < \mu < \frac{1}{3}$ | $[0, \frac{\pi}{6} - \frac{\varphi}{2}] \cup [\frac{\pi}{6} + \frac{\varphi}{2}, \frac{\pi}{2} - \frac{\varphi}{2}]$ |
| $\frac{1}{3} \leq \mu \leq 3$ | $[0, \frac{\pi}{2}]$ (all orientations) |
| $\mu > 3$ | $[\frac{\varphi}{2}, \frac{\pi}{3} - \frac{\varphi}{2}] \cup [\frac{\pi}{3} + \frac{\varphi}{2}, \frac{\pi}{2}]$ |

where $\nu(\mu) = \frac{\mu+1}{2|\mu-1|}$ and $\varphi(\mu) = \arccos(\nu(\mu))$.

The peak positions satisfy $\operatorname{span}\{y_1, \ldots, y_M\} \subseteq U$, the unstable subspace of $H_V(x_0)$.

### 2. Perturbation System ($R > 0$)

- **Explicit solution at $R = 0$ (Proposition 4.1).** For the equilateral triangle at $\theta = 0$, define $\beta = \frac{1}{2}\ln(\alpha_{12}/\alpha_{23})$ and $\gamma = 2\beta/(3\sqrt{3})$. Then:

$$\delta y_1^0 = -\gamma\, \mathbf{r}_2, \qquad \delta y_{2,3}^0 = \mp\beta\, \mathbf{r}_1 + \frac{\gamma}{2}\,\mathbf{r}_2$$

solves $F_i(\delta y, 0) = 0$ and satisfies the center-of-mass condition $\sum_i \delta y_i^0 = 0$.

- **Four-dimensional kernel (Proposition 4.2).** The linearization $D_{\delta y} F_i(\delta y^0, 0)$ has kernel $\Omega = \{(w_1, w_2, w_3) \in (\mathbb{R}^2)^3 : \langle w_k - w_i, y_k - y_i \rangle = 0\}$, decomposing as:
  - 2 translational modes $\{(v,v,v) : v \in \mathbb{R}^2\}$
  - 1 rotational mode $w_{\mathrm{rot}}$
  - 1 breathing mode $w_{\mathrm{br}} = (y_1^\perp, y_2^\perp, y_3^\perp)$

### 3. Kernel Resolution and Gauge-Fixing

- **Projection decomposition.** The system is split into perpendicular ($P_\perp F$) and parallel ($P_w F$, $w \in \Omega$) projections.
- **Gauge-fixing constraints:**
  - **(G1)** $\delta y_1 \perp \mathbf{r}_1$ — first peak has no displacement along $\mathbf{r}_1$.
  - **(G2)** $(\delta y_3 - \delta y_2) \parallel \mathbf{r}_1$ — base of triangle displaces only along $\mathbf{r}_1$.
- **Rotational equivariance (Proposition 5.1).** The one-parameter family $y_i(\theta) = R(\theta) y_i'$ generates a smooth family $\delta y^0(\theta)$ of solutions at $R = 0$.
- **Center-of-mass identity (Proposition 5.2).** $\sum_i F_i(\delta y, R) = H_V \frac{\sum_i \delta y_i}{\tilde{m}_R}$, so $\sum_i F_i = 0$ iff $\sum_i \delta y_i = 0$.

### Open Problems

1. **Prove Conjecture 5.1 (gauge-fixed invertibility).** Compute the $4 \times 4$ linearized operator on the gauge-fixed subspace explicitly and verify nonsingularity for $\mu > 3$.
2. **Extend to all admissible $\theta$.** The perturbation analysis is at $\theta = 0$; extending to $R > 0$ for all admissible $\theta$ simultaneously requires uniform estimates.
3. **Stability.** Multi-peak branches with $M \geq 2$ have at least two negative directions, implying orbital instability. Adapting asymptotic stability techniques to multi-peak states remains open.

## Repository Structure

```
nls-research/
├── report/                              # Main research paper
│   ├── Research.lyx                     # Primary LyX source
│   ├── Research.tex                     # Exported LaTeX source (1400 lines)
│   ├── Research.pdf                     # Compiled PDF
│   ├── Research.bib                     # BibTeX bibliography
│   └── archive/                         # Archived drafts
│
├── figures/                             # Figure assets (PNG) and generation
│   ├── figure_1_admissibility_binding.png
│   ├── figure_2_peak_configurations.png
│   ├── figure_3_alpha_vs_theta_vs_mu.png
│   ├── figure_4_profile_decomposition.png
│   ├── figure_5_concentration_vs_peaks.png
│   ├── figure_6_interaction_decay.png
│   ├── figure_7_methodology_flowchart.png
│   └── figures.py                       # Figure generation script
│
├── triangle/                            # Mathematica computations
│   ├── Research-Triangles.nb            # Triangle configuration analysis
│   ├── Research-Triangle-Pertubation.nb # Perturbation analysis
│   ├── Research-Triangle-Pertubation.pdf
│   └── triangle_soln.md                 # Markdown derivation of admissibility
│
├── dev/                                  # Development and verification
│   ├── notes_from_claude.tex            # Development notes
│   └── audit/                            # Mathematica audit notebooks
│       ├── NLS_Mathematical_Audit.nb
│       ├── NLS_Mathematical_Audit.m
│       └── NLS_Mathematical_Audit (2).nb
│
├── requirements.txt                     # Python dependencies
├── .forgejo/workflows/                  # CI: GitHub mirror
├── AGENTS.md                            # Agent configuration
└── .gitignore
```

## Figures

All 7 figures are in `figures/` as PNG files.

| # | Filename | Paper Reference | Description |
|---|----------|----------------|-------------|
| 1 | `figure_1_admissibility_binding` | Fig. 6 | $(\mu, \theta)$ admissibility phase diagram, binding constraints, and $\theta$-ranges vs $\mu$ |
| 2 | `figure_2_peak_configurations` | Fig. 3 | Schematic configurations: collinear ($M=2$, $M=5$), isosceles, equilateral, and perturbation $\delta y^0$ |
| 3 | `figure_3_alpha_vs_theta_vs_mu` | Fig. 5 | Interaction coefficients $\alpha_{12}$, $\alpha_{13}$, $\alpha_{23}$ vs $\theta$ for $\mu \in \{1/5, 1/2, 1, 2, 3, 5\}$ |
| 4 | `figure_4_profile_decomposition` | Fig. 2 | Concentration-compactness decomposition for $M=3$ peaks with remainder $h_R$ |
| 5 | `figure_5_concentration_vs_peaks` | Fig. 1 | Evolution of peak structure with increasing $E$; $Rz_i \to x_0$ convergence |
| 6 | `figure_6_interaction_decay` | Fig. 4 | Semilog plot of $Q_R(d)$ decay and normalized $\alpha_{ki}$ vs distance ratio |
| 7 | `figure_7_methodology_flowchart` | Fig. 1 | Analytical pipeline: Kirr-Natarajan (2018) framework → Kirr-Knox (2026) contributions |

## Paper Outline

The paper (`report/Research.tex`) is organized as follows:

1. **Introduction** — Physical motivation from optical fiber physics (NLS, Kerr nonlinearity, refractive-index profiles as external potentials); overview of contributions.
2. **Setup and Preliminaries**
   - 2.1 Hypotheses (V1)–(V4) on the potential
   - 2.2 Rescaling and concentration compactness decomposition
   - 2.3 Lyapunov-Schmidt reduction (two lemmas)
   - 2.4 The reduced algebraic system; interaction functional $Q_R$; coefficients $\alpha_{ki}$; interaction strength $\chi$; corollary on unstable subspace
3. **Classification of Limiting Configurations** — Theorem 3.1 with proofs for collinear, isosceles, equilateral, and rotational equilateral cases; admissibility analysis with unified constraint; Table 1 summary.
4. **The Perturbation System**
   - 4.1 Formulation of $F_i(\delta y, R) = 0$
   - 4.2 Explicit solution at $R = 0$ (Proposition 4.1)
   - 4.3 Kernel of the linearization (Proposition 4.2): 4D kernel with translational, rotational, and breathing modes
5. **Kernel Resolution and Gauge-Fixing**
   - 5.1 Projection decomposition ($P_\perp F$, $P_w F$)
   - 5.2 Gauge-fixing constraints (G1)–(G2)
   - 5.3 Rotational equivariance at $R = 0$ (Proposition 5.1, Corollary 5.1)
   - 5.4 Extension to $R > 0$: Conjecture 5.1 (gauge-fixed invertibility), Proposition 5.2 (center-of-mass identity)
6. **Conclusion and Open Problems** — Three concrete open problems
7. **Appendix: Notation Reference** — Comprehensive notation tables (10 subsections)

## Prerequisites

### Python

- Python 3.8+
- numpy
- matplotlib

### LaTeX / LyX

- LyX 2.4+ (for editing `.lyx` source)
- A standard LaTeX distribution (TeX Live or MiKTeX) with packages: `extarticle`, `geometry`, `amsmath`, `amsthm`, `amssymb`, `graphicx`, `booktabs`, `varwidth`, `tabularx`, `units`, `float`, `enumitem`, `babel`

### Mathematica

- Wolfram Mathematica (for `triangle/*.nb` notebooks)

## Setup and Usage

### Generating Figures

```bash
cd figures
pip install -r ../requirements.txt
python figures.py
```

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

The propagation of intense light through optical fibers is governed by the interplay between linear dispersion and Kerr nonlinearity ($\chi^{(3)}$). Starting from Maxwell's equations in a source-free, non-magnetic dielectric, the electric field is written as a slowly varying envelope modulating a carrier:

$$\mathbf{E}(\mathbf{r},t) = \tfrac{1}{2}\hat{\mathbf{x}}\bigl[F(x,y)\,A(z,t)\,e^{i(\beta_0 z - \omega_0 t)} + \text{c.c.}\bigr]$$

The slowly varying envelope approximation yields the NLS equation:

$$i\partial_z A - \frac{\beta_2}{2}\partial_T^2 A + \gamma |A|^2 A = 0$$

Spatial variations in the refractive index (fiber Bragg gratings, graded-index fibers, photonic-crystal fibers) enter as an external potential $V(\mathbf{x})$. The transverse modal equation is mathematically identical to the stationary Schrödinger equation $-\Delta\psi + V(x)\psi = E\psi$, making guided modes correspond to bound states of the effective potential.

### Mathematical Framework

The analysis builds on:

- **Concentration compactness** (Lions, 1984) — classifying minimizing sequences into compactness, vanishing, or dichotomy on unbounded domains.
- **Global bifurcation analysis** (Kirr and Natarajan, 2018) — rescaling and multi-peak decomposition of ground states via equivariant bifurcation theory; the interaction functional $Q_R$ and its exponential decay.
- **Lyapunov-Schmidt reduction** — reduction to a finite-dimensional algebraic system governing peak positions, with projection onto the translational kernel.
- **Grillakis-Shatah-Strauss stability theory** — orbital instability of multi-peak states in Hamiltonian systems with symmetry.
- **Floer-Weinstein semiclassical analysis** (1986) — concentration near non-degenerate critical points for bounded potentials.
- **Symmetry-breaking bifurcation** (Kirr, Kevrekidis, Pelinovsky, 2011) — formal hypotheses (V1)–(V4) for symmetric potentials.

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
