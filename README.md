# Multi-Peak Ground-State Solutions for the Nonlinear Schrödinger Equation with External Potential

Mathematical research on multi-peak ground-state solutions of the stationary NLS with external potential in the high-energy regime ($E \to \infty$), using concentration compactness and Lyapunov–Schmidt reduction. The work establishes existence and classification of ground-state configurations with multiple concentration peaks at a non-degenerate local maximum of the potential $V$.

## Authors

- **Eduard Kirr** — Associate Professor, University of Illinois at Urbana-Champaign
- **Ethan Knox** — Graduate Student, University of Illinois at Urbana-Champaign

## Overview

We study the stationary NLS equation with external potential in the high-energy regime:

$$\left(-\Delta + V + E\right)\psi + \sigma\left|\psi\right|^{2p}\psi = 0$$

as $E \to \infty$. The work is motivated by optical fiber physics — the envelope equation in silica fibers with Kerr nonlinearity ($\chi^{\left(3\right)}$) maps to this stationary Schrödinger equation, where the refractive-index profile becomes an external potential $V\left(\mathbf{x}\right)$.

### Hypotheses on the Potential

Following Kirr, Kevrekidis, and Pelinovsky (2011), Floer and Weinstein (1986), and Kirr and Natarajan (2018):

- **(V1)** $V \in L^\infty\left(\mathbb{R}^n\right)$ — bounded and measurable.
- **(V2)** $\lim\_{\left|x\right|\to\infty} V\left(x\right) = 0$ — vanishes at infinity.
- **(V3)** $V$ is $C^2$ near $x\_0$, a non-degenerate critical point: $\nabla V\left(x\_0\right) = 0$ and $H\_V\left(x\_0\right)$ invertible.
- **(V4)** All eigenvalues of $H\_V\left(x\_0\right)$ are strictly negative — $x\_0$ is a strict local maximum.

### Key Framework

1. **Rescaling.** Define $u\_E\left(x\right) = E^{-\frac{1}{2p}}\psi\_E\left(E^{-\frac12}x\right)$, which transforms the problem into a family of equations parametrized by $R = E^{-\frac12} \to 0$.
2. **Concentration compactness decomposition.** Each solution decomposes into $M$ well-separated profiles $u\_R$ centered at distinct points $z\_1\left(R\right), \ldots, z\_M\left(R\right)$, with a remainder $h\_R$:

$$u\_E\left(x\right) = \sum\_{i=1}^{M} u\_R\left(x - z\_i\left(R\right)\right) + h\_R\left(x\right)$$

3. **Lyapunov–Schmidt reduction.** Projection onto the kernel of the linearized operator yields a finite-dimensional algebraic system for peak positions:

$$H\_V\left(x\_0\right) y\_i = \chi \sum\_{k \in N\_i} \alpha\_{ki}\left(y\_k - y\_i\right)$$

where $H\_V\left(x\_0\right)$ is the Hessian of $V$ at the critical point $x\_0$, $\alpha\_{ki} = \frac{Q\_R\left(\left|z\_k - z\_i\right|\right)}{Q\_R\left(m\_R\right)} \in \left[0,1\right]$ are normalized interaction coefficients, and $\chi > 0$ is the overall interaction strength.

4. **Perturbation system.** The full problem at finite energy $R > 0$ is captured by $F\_i\left(\delta y, R\right) = 0$, with $\delta y$ denoting peak displacements from the $R = 0$ limiting configuration.

## Main Results

### 1. Classification of Limiting Configurations ($R = 0$)

**Theorem 3.1.** Under Hypotheses (V1)–(V4), the admissible solutions of the algebraic system with all $\alpha\_{ki} \geq 0$ are:

- **(a) Collinear** ($M \geq 2$) — Peaks aligned along a single eigenvector $\mathbf{r}\_1$ of $H\_V\left(x\_0\right)$.
- **(b) Isosceles triangle** ($M = 3$) — Peaks with $\alpha\_{12} = \alpha\_{13} = -\frac{\lambda\_2}{3}$; requires the eigenvalue constraint $\lambda\_2 = 3\lambda\_1$.
- **(c) Equilateral triangle** ($M = 3$) — All peaks at unit mutual distance; $\alpha\_{12} = \alpha\_{13} = -\frac{\lambda\_2}{3}$, $\alpha\_{23} = \frac{\lambda\_2}{6} - \frac{\lambda\_1}{2}$; requires $\frac{\left|\lambda\_2\right|}{\left|\lambda\_1\right|} \geq 3$.
- **(d) Rotational equilateral** ($M = 3$, angle $\theta$) — Equilateral triangle rotated by $\theta$ relative to the eigenvector basis, with admissible $(\mu, \theta)$ regions given below.

**Admissibility regions** (where $\mu = \frac{\left|\lambda\_2\right|}{\left|\lambda\_1\right|}$):

| $\mu$-Range | $\theta$-Range |
| --- | --- |
| $0 < \mu < \frac{1}{3}$ | $\left[0, \frac{\pi}{6} - \frac{\varphi}{2}\right] \cup \left[\frac{\pi}{6} + \frac{\varphi}{2}, \frac{\pi}{2} - \frac{\varphi}{2}\right]$ |
| $\frac{1}{3} \leq \mu \leq 3$ | $\left[0, \frac{\pi}{2}\right]$ (all orientations) |
| $\mu > 3$ | $\left[\frac{\varphi}{2}, \frac{\pi}{3} - \frac{\varphi}{2}\right] \cup \left[\frac{\pi}{3} + \frac{\varphi}{2}, \frac{\pi}{2}\right]$ |

where $\nu\left(\mu\right) = \frac{\mu+1}{2\left|\mu-1\right|}$ and $\varphi\left(\mu\right) = \arccos\left(\nu\left(\mu\right)\right)$.

The peak positions satisfy $\text{span}\lbrace y\_1, \ldots, y\_M\rbrace \subseteq U$, the unstable subspace of $H\_V\left(x\_0\right)$.

### 2. Perturbation System ($R > 0$)

- **Explicit solution at $R = 0$ (Proposition 4.1).** For the equilateral triangle at $\theta = 0$, define $\beta = \frac{1}{2}\ln\left(\frac{\alpha\_{12}}{\alpha\_{23}}\right)$ and $\gamma = \frac{2\beta}{3\sqrt{3}}$. Then:

$$\delta y\_1^0 = -\gamma\, \mathbf{r}\_2, \qquad \delta y\_{2,3}^0 = \mp\beta\, \mathbf{r}\_1 + \frac{\gamma}{2}\,\mathbf{r}\_2$$

solves $F\_i\left(\delta y, 0\right) = 0$ and satisfies the center-of-mass condition $\sum\_i \delta y\_i^0 = 0$.

- **Four-dimensional kernel (Proposition 4.2).** The linearization $D\_{\delta y} F\_i\left(\delta y^0, 0\right)$ has kernel $\Omega = \lbrace\left(w\_1, w\_2, w\_3\right) \in \left(\mathbb{R}^2\right)^3 : \left \langle w\_k - w\_i, y\_k - y\_i \right \rangle = 0\rbrace$, decomposing as:
  - 2 translational modes $\lbrace\left(v,v,v\right) : v \in \mathbb{R}^2\rbrace$
  - 1 rotational mode $w\_{\mathrm{rot}}$
  - 1 breathing mode $w\_{\mathrm{br}} = \left(y\_1^\perp, y\_2^\perp, y\_3^\perp\right)$

### 3. Kernel Resolution and Gauge-Fixing

- **Projection decomposition.** The system is split into perpendicular ($P\_\perp F$) and parallel ($P\_w F$, $w \in \Omega$) projections.
- **Gauge-fixing constraints:**
  - **(G1)** $\delta y\_1 \perp \mathbf{r}\_1$ — first peak has no displacement along $\mathbf{r}\_1$.
  - **(G2)** $\left(\delta y\_3 - \delta y\_2\right) \parallel \mathbf{r}\_1$ — base of triangle displaces only along $\mathbf{r}\_1$.
- **Rotational equivariance (Proposition 5.1).** The one-parameter family $y\_i\left(\theta\right) = R\left(\theta\right) y\_i'$ generates a smooth family $\delta y^0\left(\theta\right)$ of solutions at $R = 0$.
- **Center-of-mass identity (Proposition 5.2).** $\sum\_i F\_i\left(\delta y, R\right) = H\_V \frac{\sum\_i \delta y\_i}{\tilde{m}\_R}$, so $\sum\_i F\_i = 0$ iff $\sum\_i \delta y\_i = 0$.

### Open Problems

1. **Prove Conjecture 5.1 (gauge-fixed invertibility).** Compute the $4 \times 4$ linearized operator on the gauge-fixed subspace explicitly and verify nonsingularity for $\mu > 3$.
2. **Extend to all admissible $\theta$.** The perturbation analysis is at $\theta = 0$; extending to $R > 0$ for all admissible $\theta$ simultaneously requires uniform estimates.
3. **Stability.** Multi-peak branches with $M \geq 2$ have at least two negative directions, implying orbital instability. Adapting asymptotic stability techniques to multi-peak states remains open.

## Repository Structure

```
nls-research/
├── report/                                  # Main research paper
│   ├── Research.lyx                         # Primary LyX source (source of truth, ~7400 lines)
│   ├── Research.tex                         # LyX-exported LaTeX (~1200 lines)
│   ├── Research.pdf                         # Compiled PDF
│   └── Research.bib                         # BibTeX bibliography (14 entries)
│
├── figures/                                 # Figure generation and output
│   ├── figures.py                           # Single-file script generating all figures (~1400 lines)
│   ├── figure_1_admissibility_binding.png
│   ├── figure_2_peak_configurations.png
│   ├── figure_2b_perturbation.png
│   ├── figure_3_alpha_vs_theta_vs_mu.png
│   ├── figure_4_profile_decomposition.png
│   ├── figure_5_concentration_vs_peaks.png
│   ├── figure_6_interaction_decay.png
│   └── figure_7_methodology_flowchart.png
│
├── requirements.txt                         # Python dependencies (numpy, matplotlib)
├── .forgejo/workflows/                      # CI: Forgejo → GitHub mirror
└── .gitignore
```

## Figures

All figures are generated by `figures/figures.py` and output as PNG at 300 DPI. The `create_figure_2` function produces both Figure 2 and Figure 2b.

| # | Filename | Description |
| --- | --- | --- |
| 1 | `figure_1_admissibility_binding` | $(\mu, \theta)$ admissibility phase diagram with binding constraints |
| 2 | `figure_2_peak_configurations` | Schematic configurations: collinear ($M=2$, $M=5$), isosceles, and equilateral triangles |
| 2b | `figure_2b_perturbation` | Perturbation displacement $\delta y^0$ for the equilateral triangle at $\theta = 0$ |
| 3 | `figure_3_alpha_vs_theta_vs_mu` | Interaction coefficients $\alpha\_{12}$, $\alpha\_{13}$, $\alpha\_{23}$ vs $\theta$ for various $\mu$ |
| 4 | `figure_4_profile_decomposition` | Concentration-compactness decomposition for $M=3$ peaks with remainder $h\_R$ |
| 5 | `figure_5_concentration_vs_peaks` | Evolution of peak structure with increasing $E$; $Rz\_i \to x\_0$ convergence |
| 6 | `figure_6_interaction_decay` | Semilog plot of $Q\_R(d)$ decay and normalized $\alpha\_{ki}$ vs distance ratio |
| 7 | `figure_7_methodology_flowchart` | Analytical pipeline: Kirr–Natarajan (2018) framework → Kirr–Knox (2026) contributions |

## Paper Outline

The report (`Research.lyx`, exported as `Research.tex`) uses the `extarticle` document class and is organized as follows:

1. **Introduction** — Physical motivation from optical fiber physics (NLS, Kerr nonlinearity, refractive-index profiles as external potentials); overview of contributions.
2. **Setup and Preliminaries**
   - 2.1 Hypotheses (V1)–(V4) on the potential
   - 2.2 Rescaling and concentration compactness decomposition
   - 2.3 Lyapunov–Schmidt reduction (two lemmas)
   - 2.4 The reduced algebraic system; interaction functional $Q\_R$; coefficients $\alpha\_{ki}$; interaction strength $\chi$; corollary on unstable subspace
3. **Classification of Limiting Configurations** — Theorem 3.1 with proofs for each configuration type
   - 3.1 Collinear configurations
   - 3.2 Isosceles triangle
   - 3.3 Equilateral triangle
   - 3.4 Rotational equilateral triangle
   - 3.5 Admissibility analysis with unified constraint; summary table
4. **The Perturbation System**
   - 4.1 Formulation of $F\_i\left(\delta y, R\right) = 0$
   - 4.2 Explicit solution at $R = 0$ (Proposition 4.1)
   - 4.3 Kernel of the linearization (Proposition 4.2): 4D kernel with translational, rotational, and breathing modes
5. **Kernel Resolution and Gauge-Fixing**
   - 5.1 Projection decomposition ($P\_\perp F$, $P\_w F$)
   - 5.2 Gauge-fixing constraints (G1)–(G2)
   - 5.3 Rotational equivariance at $R = 0$ (Proposition 5.1, Corollary 5.1)
   - 5.4 Extension to $R > 0$: Conjecture 5.1 (gauge-fixed invertibility), Proposition 5.2 (center-of-mass identity)
6. **Conclusion, Open Problems, and Acknowledgements** — Summary of contributions; three concrete open problems
7. **Notation Reference** — Comprehensive notation tables

## Prerequisites

### Python

- Python 3.8+
- numpy
- matplotlib

### LaTeX / LyX

- LyX 2.4+ (for editing the `.lyx` source)
- A standard LaTeX distribution (TeX Live or MiKTeX) with packages: `extarticle`, `geometry`, `amsmath`, `amsthm`, `amssymb`, `graphicx`, `booktabs`, `tabularx`, `units`, `float`

### Mathematica

- Wolfram Mathematica (for `triangle/*.nb` notebooks)

## Setup and Usage

### Generating Figures

```bash
cd figures
pip install -r ../requirements.txt
python figures.py
```

Figures are written to `figures/` as PNG files at 300 DPI.

### Compiling the Paper

From `report/`:

```bash
pdflatex Research.tex
bibtex Research
pdflatex Research.tex
pdflatex Research.tex
```

Or open `Research.lyx` in LyX and use the built-in export/compile. Note that `Research.lyx` is the source of truth — `Research.tex` is generated from it via LyX export.

## Background and Context

### Physical Motivation

The propagation of intense light through optical fibers is governed by the interplay between linear dispersion and Kerr nonlinearity ($\chi^{(3)}$). Starting from Maxwell's equations in a source-free, non-magnetic dielectric, the electric field is written as a slowly varying envelope modulating a carrier:

$$\mathbf{E}\left(\mathbf{r},t\right) = \frac12\hat{\mathbf{x}}\left[F\left(x,y\right)\,A\left(z,t\right)\,e^{i\left(\beta\_0 z - \omega\_0 t\right)} + \text{c.c.}\right]$$

The slowly varying envelope approximation yields the NLS equation:

$$i\partial\_z A - \frac{\beta\_2}{2}\partial\_T^2 A + \gamma \left|A\right|^2 A = 0$$

Spatial variations in the refractive index (fiber Bragg gratings, graded-index fibers, photonic-crystal fibers) enter as an external potential $V\left(\mathbf{x}\right)$. The transverse modal equation is mathematically identical to the stationary Schrödinger equation $-\Delta\psi + V\left(x\right)\psi = E\psi$, making guided modes correspond to bound states of the effective potential.

### Mathematical Framework

The analysis builds on:

- **Concentration compactness** (Lions, 1984) — classifying minimizing sequences into compactness, vanishing, or dichotomy on unbounded domains.
- **Global bifurcation analysis** (Kirr and Natarajan, 2018) — rescaling and multi-peak decomposition of ground states via equivariant bifurcation theory; the interaction functional $Q\_R$ and its exponential decay.
- **Lyapunov–Schmidt reduction** — reduction to a finite-dimensional algebraic system governing peak positions, with projection onto the translational kernel.
- **Grillakis–Shatah–Strauss stability theory** — orbital instability of multi-peak states in Hamiltonian systems with symmetry.
- **Floer–Weinstein semiclassical analysis** (1986) — concentration near non-degenerate critical points for bounded potentials.
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