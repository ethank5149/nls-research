"""
Figure 7 — Analytical Pipeline: PDE → Multi-Peak Ground States
Compact two-column layout for inclusion in extarticle report.
"""
import warnings
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'font.size': 7.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

C_BLUE_FG  = '#2c5f8a';  C_BLUE_BG  = '#dce8f2'
C_GRAY_FG  = '#4a4a4a';  C_GRAY_BG  = '#edeceb'
C_GREEN_FG = '#2a6b4a';  C_GREEN_BG = '#d6ead6'
C_RED_FG   = '#9b2d2d';  C_RED_BG   = '#f5e0e0'
C_OPEN_FG  = '#7a6520';  C_OPEN_BG  = '#faf3dc'
C_ARROW    = '#3d3d3d'
C_HDR_R    = '#7b3333'


def box(ax, x, y, w, h, text, fc, ec, fs=7, bold=False, lw=1.2, ls='-'):
    b = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle="round,pad=0.08",
        facecolor=fc, edgecolor=ec, lw=lw, linestyle=ls,
    )
    ax.add_patch(b)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs,
            fontweight='bold' if bold else 'normal',
            multialignment='center', color='#1a1a1a', linespacing=1.25)


def arr(ax, x1, y1, x2, y2, color=C_ARROW, lw=1.0, ls='-'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle=ls, shrinkA=1, shrinkB=1))


def col_header(ax, x, y, w, text, color):
    rect = FancyBboxPatch(
        (x - w/2, y - 0.22), w, 0.44,
        boxstyle="round,pad=0.06",
        facecolor=color, edgecolor=color, lw=1.5, alpha=0.90,
    )
    ax.add_patch(rect)
    ax.text(x, y, text, ha='center', va='center',
            fontsize=7.5, fontweight='bold', color='white')


def create():
    print("Creating Figure 7 \u2013 Methodology Flowchart...", end=' ', flush=True)

    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    ax.set_xlim(-0.2, 12.2)
    ax.set_ylim(1.5, 17.5)
    ax.axis('off')
    ax.set_aspect('equal')

    Lx  = 2.75;  Rx  = 9.25
    Bw  = 4.6
    sp  = 1.12;  bh  = 0.62;  bhL = 0.74

    Y_top = 17.0

    col_header(ax, Lx, Y_top, Bw + 0.2, 'Kirr\u2013Natarajan (2018)', C_BLUE_FG)
    col_header(ax, Rx, Y_top, Bw + 0.2, 'Kirr\u2013Knox (2026)', C_HDR_R)

    # ═══════════════ LEFT COLUMN ════════════════════════════════════
    Y = Y_top - 0.80

    box(ax, Lx, Y, Bw, bh,
        r'$(-\Delta + V + E)\psi + \sigma|\psi|^{2p}\psi = 0$'
        '\nFull PDE on ' r'$\mathbb{R}^n$',
        C_BLUE_BG, C_BLUE_FG, fs=6.8, bold=True)

    Y -= sp; arr(ax, Lx, Y+sp-bh/2, Lx, Y+bh/2)
    box(ax, Lx, Y, Bw, bh,
        r'Rescaling:  $u_E = E^{-1/(2p)}\psi_E(E^{-1/2}\cdot)$'
        '\n' r'$\|u_E\|_{H^1}$ bounded;  $R = E^{-1/2}\to 0$',
        C_GRAY_BG, C_GRAY_FG, fs=6.3)

    Y -= sp; arr(ax, Lx, Y+sp-bh/2, Lx, Y+bh/2)
    box(ax, Lx, Y, Bw, bh,
        'Conc.-Compactness (Lions)'
        '\n' r'$u_E \approx \sum_{i=1}^{M} u_R(\cdot - z_i) + h_R$',
        C_GRAY_BG, C_GRAY_FG, fs=6.3)

    Y -= sp; arr(ax, Lx, Y+sp-bh/2, Lx, Y+bh/2)
    box(ax, Lx, Y, Bw, bh,
        'Lyapunov\u2013Schmidt (Lemma 2)'
        '\n' r'$P_\perp F = 0 \Rightarrow h_R = h_R(z_1,\ldots,z_M,R)$',
        C_GRAY_BG, C_GRAY_FG, fs=6.3)

    Y -= sp*1.05; arr(ax, Lx, Y+sp*1.05-bh/2, Lx, Y+bhL/2)
    box(ax, Lx, Y, Bw, bhL,
        'Reduced System'
        '\n' r'$H_V y_i = \sum_{k \in N_i} \alpha_{ki}\,\chi(y_k - y_i)$'
        '\n' r'peaks $\in$ unstable manifold of $H_V$',
        C_GREEN_BG, C_GREEN_FG, fs=6.3, bold=True)

    Y_fork = Y - bhL/2 - 0.12

    Y_res = Y_fork - 0.55
    arr(ax, Lx, Y_fork, Lx, Y_res + bhL/2)
    box(ax, Lx, Y_res, Bw, bhL,
        'Thm 4.3: one profile per critical point'
        '\n' r'$\|U_E - \sum \mu_i u_\infty(\cdot - x_i\sqrt{E})\|_{H^2} < \varepsilon$'
        '\nExistence + uniqueness + stability',
        C_GREEN_BG, C_GREEN_FG, fs=6.3)

    Y_tag = Y_res - bhL/2 - 0.35
    arr(ax, Lx, Y_res - bhL/2, Lx, Y_tag + 0.14)
    box(ax, Lx, Y_tag, Bw*0.75, 0.28,
        'Resolved for distinct critical points',
        C_GREEN_BG, '#1a5632', fs=6, bold=True)

    Y_open = Y_tag - 0.14 - 0.55
    arr(ax, Lx, Y_tag - 0.14, Lx, Y_open + bhL/2)
    box(ax, Lx, Y_open, Bw, bhL,
        'Open Problem (K-N Remark 4.1):'
        '\n' r'Multiple profiles $\to$ same critical point'
        '\n' r'along eigenvectors with $\lambda_j < 0$',
        C_RED_BG, C_RED_FG, fs=6.3, bold=True, lw=1.6)

    # Left border
    pad = 0.18
    yb_L = Y_open - bhL/2 - 0.10
    yt_L = Y_top - 0.48
    left_rect = FancyBboxPatch(
        (Lx - Bw/2 - pad, yb_L), Bw + 2*pad, yt_L - yb_L,
        boxstyle="round,pad=0.10",
        facecolor='none', edgecolor=C_BLUE_FG, lw=0.7, alpha=0.22,
    )
    ax.add_patch(left_rect)

    # ═══════════════ BRIDGE ═════════════════════════════════════════
    Y_H = Y_top - 0.80
    arr(ax, Lx + Bw/2 + 0.08, Y_open + 0.05,
        Rx - Bw/2 - 0.08, Y_H - 0.05,
        color=C_RED_FG, lw=1.5)
    mid_x = (Lx + Bw/2 + Rx - Bw/2) / 2
    mid_y = (Y_open + Y_H) / 2
    ax.text(mid_x + 0.1, mid_y + 0.15, 'this work',
            ha='center', va='bottom', fontsize=5.5, fontstyle='italic',
            color=C_RED_FG)

    # ═══════════════ RIGHT COLUMN ═══════════════════════════════════
    Y = Y_H

    box(ax, Rx, Y, Bw, bh,
        r'Classify Limiting Configs ($R = 0$)'
        '\nCollinear | Isosceles | Equilateral | Rot. ' r'$\theta$',
        C_GRAY_BG, C_GRAY_FG, fs=6.3)

    Y -= sp; arr(ax, Rx, Y+sp-bh/2, Rx, Y+bh/2)
    box(ax, Rx, Y, Bw, bh,
        r'Admissibility:  $\alpha_{ki} \geq 0$'
        '\n' r'$\frac{1}{3} \leq \mu \leq 3$: all $\theta$;  restricted otherwise',
        C_GREEN_BG, C_GREEN_FG, fs=6.3)

    Y -= sp; arr(ax, Rx, Y+sp-bh/2, Rx, Y+bh/2)
    box(ax, Rx, Y, Bw, bh,
        r'Perturbation System:  $F_i(\delta y, R) = 0$'
        '\n' r'Extended to $R = 0$ by continuity',
        C_BLUE_BG, C_BLUE_FG, fs=6.3)

    Y -= sp*1.05; arr(ax, Rx, Y+sp*1.05-bh/2, Rx, Y+bhL/2)
    box(ax, Rx, Y, Bw, bhL,
        r'Explicit Solution at $R=0$ (Claim 3)'
        '\n' r'$\delta y_1^0=-\gamma\mathbf{r}_2$,  '
        r'$\delta y_{2,3}^0=\mp\beta\mathbf{r}_1+\frac{\gamma}{2}\mathbf{r}_2$',
        C_BLUE_BG, C_BLUE_FG, fs=6.3)

    Y -= sp*1.05; arr(ax, Rx, Y+sp*1.05-bhL/2, Rx, Y+bhL/2)
    box(ax, Rx, Y, Bw, bhL,
        'Kernel Obstruction (Claim 4)'
        '\n' r'$\Omega = \ker(D_{\delta y}F_i) \supseteq \{(v,v,v):v \in \mathbb{R}^2\}$'
        '\nIFT blocked by translational invariance',
        C_RED_BG, C_RED_FG, fs=6.3, bold=True, lw=1.6)

    Y -= sp*1.05; arr(ax, Rx, Y+sp*1.05-bhL/2, Rx, Y+bhL/2)
    box(ax, Rx, Y, Bw, bhL,
        r'Decompose:  $\mathbb{R}^{2\times 3}=\Omega\oplus\Omega^\perp$'
        '\n' r'$P_\perp F$: IFT on $\Omega^\perp$  |  '
        r'$P_w F$: via equivariance',
        C_GRAY_BG, C_GRAY_FG, fs=6.0)

    Y -= sp; arr(ax, Rx, Y+sp-bhL/2, Rx, Y+bh/2)
    box(ax, Rx, Y, Bw, bh,
        r'Gauge-Fixing:  $\delta y_1 \perp \mathbf{r}_1$  ;  '
        r'$(\delta y_3-\delta y_2) \parallel \mathbf{r}_1$',
        C_GRAY_BG, C_GRAY_FG, fs=6.3)

    Y -= sp*1.05; arr(ax, Rx, Y+sp*1.05-bh/2, Rx, Y+bhL/2)
    box(ax, Rx, Y, Bw, bhL,
        'IFT + Lift via Lemma 2'
        '\n' r'$\Rightarrow$ exact multi-peak $\psi_E$ for $E \gg 1$'
        '\n(contingent on completing open steps)',
        C_GREEN_BG, '#1a5632', fs=6.3, bold=True, lw=1.3, ls='--')
    Y_Q = Y

    # Right border
    yb_R = Y_Q - bhL/2 - 0.10
    yt_R = Y_top - 0.48
    right_rect = FancyBboxPatch(
        (Rx - Bw/2 - pad, yb_R), Bw + 2*pad, yt_R - yb_R,
        boxstyle="round,pad=0.10",
        facecolor='none', edgecolor=C_HDR_R, lw=0.7, alpha=0.22,
    )
    ax.add_patch(right_rect)

    # Open steps callout
    Y_os = yb_R - 0.42
    box(ax, Rx, Y_os, Bw, 0.55,
        'Remaining:  (1) Verify equivariance of ' r'$F_i$'
        '\n(2) Invertibility on gauge-fixed subspace'
        '\n(3) Extend to continuous ' r'$\theta$',
        C_OPEN_BG, C_OPEN_FG, fs=5.8, lw=1.0, ls='--')
    arr(ax, Rx, Y_os + 0.28, Rx, yb_R,
        color=C_OPEN_FG, lw=0.7, ls='--')

    # Legend
    leg_y = min(yb_L, Y_os - 0.28) - 0.50
    leg_x = 0.1; dx = 3.05
    for i, (fc, ec, label) in enumerate([
        (C_BLUE_BG, C_BLUE_FG, 'Equations / systems'),
        (C_GRAY_BG, C_GRAY_FG, 'Analytical machinery'),
        (C_GREEN_BG, C_GREEN_FG, 'Results / resolved'),
        (C_RED_BG, C_RED_FG, 'Obstruction / open'),
    ]):
        xx = leg_x + i * dx
        swatch = FancyBboxPatch(
            (xx, leg_y - 0.09), 0.26, 0.18,
            boxstyle="round,pad=0.03",
            facecolor=fc, edgecolor=ec, lw=0.7,
        )
        ax.add_patch(swatch)
        ax.text(xx + 0.34, leg_y, label, fontsize=5.5,
                va='center', color='#333')

    # Adjust ylim to crop tightly to content
    ax.set_ylim(leg_y - 0.45, 17.5)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig('figures/output/figure_7_methodology_flowchart.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig('figures/output/figure_7_methodology_flowchart.png', bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()