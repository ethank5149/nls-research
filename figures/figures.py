import sys, os, warnings, matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from pathlib import Path


sys.path.insert(0, os.path.dirname(__file__))
warnings.filterwarnings("ignore")
matplotlib.use("Agg")

C_BLUE_FG = "#2c5f8a"
C_BLUE_BG = "#dce8f2"
C_GRAY_FG = "#4a4a4a"
C_GRAY_BG = "#edeceb"
C_GREEN_FG = "#2a6b4a"
C_GREEN_BG = "#d6ead6"
C_RED_FG = "#9b2d2d"
C_RED_BG = "#f5e0e0"
C_OPEN_FG = "#7a6520"
C_OPEN_BG = "#faf3dc"
C_ARROW = "#3d3d3d"
C_HDR_R = "#7b3333"

OUTPUT_DIR = Path(__file__).resolve().parent


# ─────────── Standardized rcParams for publication-quality figures ───────────
RC_PARAMS = {
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
}


def apply_rcparams():
    """Apply the standardized rcParams."""
    import warnings

    warnings.filterwarnings("ignore")
    plt.rcParams.update(RC_PARAMS)


apply_rcparams()
# ─────────────────────────────────────────────────────────────────────────────


def compute_alphas(mu, theta):
    """Compute interaction coefficients α₁₂, α₁₃, α₂₃ for given μ and θ.

    Uses λ₁ = -1, λ₂ = -μ convention.
    """
    lam1, lam2 = -1.0, -mu
    s, c = np.sin(theta), np.cos(theta)
    a12 = (1 / 3) * (
        lam2 * (-np.sqrt(3) * s - c) * c + lam1 * (-s + np.sqrt(3) * c) * s
    )
    a13 = (1 / 3) * (lam2 * (np.sqrt(3) * s - c) * c - lam1 * (s + np.sqrt(3) * c) * s)
    a23 = (1 / 6) * (
        lam2 * (2 * np.cos(2 * theta) - 1) - lam1 * (2 * np.cos(2 * theta) + 1)
    )
    return a12, a13, a23


def box(ax, x, y, w, h, text, fc, ec, fs=9, bold=False, lw=1.5, ls="-"):
    """Draw a rounded box with centered text."""
    b = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.15",
        facecolor=fc,
        edgecolor=ec,
        lw=lw,
        linestyle=ls,
    )
    ax.add_patch(b)
    ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        fontweight="bold" if bold else "normal",
        multialignment="center",
    )


def arr(ax, x1, y1, x2, y2, color="#2c3e50", lw=1.5, ls="-"):
    """Draw an arrow annotation."""
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->", color=color, lw=lw, linestyle=ls, shrinkA=1, shrinkB=1
        ),
    )


def draw_peak(ax, pos, label, color="#2c3e50", ms=9, offset=(0.06, 0.10)):
    ax.plot(*pos, "o", color=color, ms=ms, zorder=5, clip_on=False)
    ax.annotate(
        label, xy=pos, xytext=(pos[0] + offset[0], pos[1] + offset[1]),
        fontsize=10, color=color, ha="center", clip_on=False,
    )


def draw_edge(ax, p1, p2, alpha_val=None, color="#7f8c8d", lw=1.3,
              label_offset_scale=0.10, label_fs=7):
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color=color, lw=lw, zorder=2)
    if alpha_val is not None:
        mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = np.sqrt(dx**2 + dy**2)
        nx, ny = -dy / length * label_offset_scale, dx / length * label_offset_scale
        ax.text(
            mx + nx, my + ny, alpha_val, fontsize=label_fs,
            ha="center", va="center", color="#2980b9", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.10", facecolor="white",
                      edgecolor="#ddd", alpha=0.9),
            zorder=8,
        )


def draw_eigvec(ax, origin, direction, label, color="#c0392b", fs=9):
    ax.annotate(
        "", xy=(origin[0] + direction[0], origin[1] + direction[1]), xytext=origin,
        arrowprops=dict(arrowstyle="->", color=color, lw=1.3),
    )
    ax.text(
        origin[0] + direction[0] * 1.18, origin[1] + direction[1] * 1.18,
        label, fontsize=fs, color=color, ha="center", va="center",
        zorder=10,
        bbox=dict(boxstyle="round,pad=0.05", facecolor="white",
                  edgecolor="none", alpha=0.85),
    )

def info_box(ax, x, y, text, fs=8, fc="#fef9e7", ec="#ddd"):
    ax.text(
        x, y, text, fontsize=fs, ha="center", va="top", color="#444",
        bbox=dict(boxstyle="round,pad=0.25", facecolor=fc, edgecolor=ec),
    )

def save_figure(fig, output_dir, stem, pad_inches=0.15, *args, **kwargs):
    """Save figure in PNG format."""
    fig.savefig(
        output_dir / f"{stem}.png",
        bbox_inches="tight",
        pad_inches=pad_inches,
        *args,
        **kwargs,
    )
    plt.close(fig)


def create_figure_1():
    """Figure 1 — Admissibility & Binding Constraint (μ,θ) Phase Diagram."""
    print(
        "Creating Figure 1 - Admissibility and Binding Constraints...",
        end=" ",
        flush=True,
    )

    mu_vals = np.linspace(0.01, 6, 1000)
    theta_vals = np.linspace(0, np.pi / 2, 1000)
    MU, TH = np.meshgrid(mu_vals, theta_vals)

    a12, a13, a23 = compute_alphas(MU, TH)
    admissible = (a12 >= -1e-12) & (a13 >= -1e-12) & (a23 >= -1e-12)

    stack = np.stack([a12, a13, a23], axis=-1)
    binding = np.argmin(stack, axis=-1).astype(float)
    binding_masked = np.where(admissible, binding, 3.0)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Refined palette: muted teal / warm gold / cool slate / light grey
    c_a12 = "#3a7ca5"   # steel blue
    c_a13 = "#d4915e"   # warm terracotta
    c_a23 = "#6b8e7b"   # sage green
    c_inad = "#e8e4df"  # warm light grey
    cmap = ListedColormap([c_a12, c_a13, c_a23, c_inad])
    ax.pcolormesh(
        MU, TH, binding_masked, cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto"
    )

    ax.contour(
        MU, TH, admissible.astype(float), levels=[0.5], colors="#1a1a2e", linewidths=2.0
    )

    ax.axvline(x=1 / 3, color="white", ls="--", lw=1.3, alpha=0.9)
    ax.axvline(x=3, color="white", ls="--", lw=1.3, alpha=0.9)
    ax.axvline(x=1, color="white", ls=":", lw=1.2, alpha=0.6)

    # "All θ admissible" label placed just below x-axis
    ax.text(
        1.65,
        -0.07,
        r"All $\theta$ admissible",
        fontsize=8.5,
        ha="center",
        va="top",
        color="#333",
        fontweight="bold",
        style="italic",
        clip_on=False,
    )
    # Bracket lines from 1/3 to 3 just below axis
    brac_y = -0.03
    ax.plot([1/3, 1/3], [0, brac_y], color="#555", lw=0.8, clip_on=False)
    ax.plot([3, 3], [0, brac_y], color="#555", lw=0.8, clip_on=False)
    ax.plot([1/3, 3], [brac_y, brac_y], color="#555", lw=0.8, clip_on=False)

    ax.plot(3, 0, "o", color="white", ms=6, mec="#1a1a2e", mew=1.5, zorder=5)
    ax.annotate(
        r"$\lambda_2 = 3\lambda_1$",
        xy=(3, 0),
        xytext=(4.2, np.pi / 6),
        fontsize=9,
        color="#1a1a2e",
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#1a1a2e", lw=1.2),
        bbox=dict(
            boxstyle="round,pad=0.2", facecolor="white", edgecolor="#ccc", alpha=0.9
        ),
    )

    ax.text(
        5.0,
        np.pi / 3,
        "Restricted\nangular\nwindows",
        fontsize=8,
        ha="center",
        color="#555",
        style="italic",
    )

    ax.set_xlabel(r"$\mu$", labelpad=6)
    ax.set_ylabel(r"$\theta$", rotation=0, labelpad=12)
    ax.set_xlim(0, 6)
    ax.set_ylim(0, np.pi / 2)
    ax.set_box_aspect(1)
    ax.set_xticks([0, 1 / 3, 3, 6])
    ax.set_xticklabels([r"$0$", r"$\frac{1}{3}$", r"$3$", r"$6$"])
    ax.set_yticks([0, np.pi / 6, np.pi / 3, np.pi / 2])
    ax.set_yticklabels(
        [r"$0$", r"$\frac{\pi}{6}$", r"$\frac{\pi}{3}$", r"$\frac{\pi}{2}$"]
    )

    legend_elements = [
        mpatches.Patch(
            color=c_a12,
            label=r"$\alpha_{12}$ binding (apex $\leftrightarrow$ base-left)",
        ),
        mpatches.Patch(
            color=c_a13,
            label=r"$\alpha_{13}$ binding (apex $\leftrightarrow$ base-right)",
        ),
        mpatches.Patch(color=c_a23, label=r"$\alpha_{23}$ binding (base edge)"),
        mpatches.Patch(
            color=c_inad, label=r"Inadmissible ($\exists\,\alpha_{ij} < 0$)"
        ),
        Line2D([0], [0], color="#1a1a2e", lw=2, label="Admissibility boundary"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=7.5,
        framealpha=0.95,
        edgecolor="#ccc",
        borderpad=0.8,
        handletextpad=0.6,
    )

    fig.tight_layout(pad=1.2)
    save_figure(fig, OUTPUT_DIR, "figure_1_admissibility_binding")
    print("Done!")


def create_figure_2():
    print("Creating Figure 2 - Peak Configurations (2x2 + separate)...", end=" ", flush=True)

    # ── Layout geometry ──
    # Manual axes placement: [left, bottom, width, height] in figure coords.
    # Two rows: top row is short (collinear), bottom row is tall (triangles).
    # All panels in a row share the exact same bottom and height.
    margin_l, margin_r = 0.03, 0.03
    margin_b, margin_t = 0.02, 0.05  # top margin for suptitle area
    col_gap = 0.04
    row_gap = 0.04

    col_w = (1.0 - margin_l - margin_r - col_gap) / 2.0
    top_h = 0.22                      # top row height (fraction of figure)
    bot_h = 1.0 - margin_b - margin_t - row_gap - top_h

    top_bottom = 1.0 - margin_t - top_h
    bot_bottom = margin_b

    left_x = margin_l
    right_x = margin_l + col_w + col_gap

    fig = plt.figure(figsize=(10, 8))

    # ── (a) Collinear M=2 ──────────────────────────────────────────────────
    ax = fig.add_axes([left_x, top_bottom, col_w, top_h])
    peaks = [(-0.5, 0), (0.5, 0)]
    draw_edge(ax, peaks[0], peaks[1], r"$\alpha_{12}$")
    draw_peak(ax, peaks[0], r"$y_{-1}$")
    draw_peak(ax, peaks[1], r"$y_1$")
    draw_eigvec(ax, (-0.15, -0.28), (0.38, 0), r"$\mathbf{r}_1$")
    ax.annotate(
        "", xy=(0.5, -0.14), xytext=(-0.5, -0.14),
        arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8),
    )
    ax.text(0, -0.21, r"$|y_1 - y_{-1}| = 1$", fontsize=8, ha="center", color="#555")
    info_box(ax, 0, -0.38,
             r"$H_V y_{\pm 1} = \alpha_{12}\mathbf{r}_1$"
             "\n"
             r"$\Rightarrow \alpha_{12} = -\lambda_1$")
    ax.set_xlim(-0.85, 0.85)
    ax.set_ylim(-0.55, 0.18)
    ax.axis("off")
    ax.set_title(r"(a) Collinear, $M=2$", fontsize=11, pad=4)

    # ── (b) Collinear M=5 ─────────────────────────────────────────────────
    ax = fig.add_axes([right_x, top_bottom, col_w, top_h])
    peaks_5 = [(-1, 0), (-0.5, 0), (0, 0), (0.5, 0), (1, 0)]
    labels_5 = [r"$y_{-2}$", r"$y_{-1}$", r"$y_0$", r"$y_1$", r"$y_2$"]
    for i in range(len(peaks_5) - 1):
        draw_edge(ax, peaks_5[i], peaks_5[i + 1])
    for p, l in zip(peaks_5, labels_5):
        draw_peak(ax, p, l)
    draw_eigvec(ax, (-0.25, -0.25), (0.38, 0), r"$\mathbf{r}_1$")
    info_box(ax, 0, -0.38,
             r"$H_V(i\mathbf{r}_1) = (\alpha_{i,i+1} + \alpha_{i,i-1})\mathbf{r}_1$"
             "\nNearest-neighbor coupling only",
             fs=7.5)
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-0.55, 0.18)
    ax.axis("off")
    ax.set_title(r"(b) Collinear, $M=5$ (odd)", fontsize=11, pad=4)

    # ── (c) Isosceles triangle ─────────────────────────────────────────────
    ax = fig.add_axes([left_x, bot_bottom, col_w, bot_h])
    a_iso, h_iso = 0.55, 0.48
    y1 = (0, h_iso)
    y2 = (-a_iso, -h_iso / 2)
    y3 = (a_iso, -h_iso / 2)
    for p1, p2, albl in [
        (y1, y2, r"$\alpha_{12}$"),
        (y1, y3, r"$\alpha_{13}$"),
        (y2, y3, r"$\alpha_{23}=0$"),
    ]:
        c_edge = "#e74c3c" if "=0" in albl else "#7f8c8d"
        draw_edge(ax, p1, p2, albl, color=c_edge, label_offset_scale=0.10)
    draw_peak(ax, y1, r"$y_1$")
    draw_peak(ax, y2, r"$y_2$", offset=(-0.12, 0.08))
    draw_peak(ax, y3, r"$y_3$")
    ev_origin = (-0.82, -0.48)
    draw_eigvec(ax, ev_origin, (0.33, 0), r"$\mathbf{r}_1$", "#c0392b")
    draw_eigvec(ax, ev_origin, (0, 0.33), r"$\mathbf{r}_2$", "#2980b9")
    ax.text(
        0.60, 0.46,
        r"$a > \frac{1}{2}$" "\n" r"$h < \frac{1}{\sqrt{3}}$",
        fontsize=8, color="#8e44ad",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="#f5eef8", edgecolor="#d2b4de"),
    )
    info_box(ax, 0, -0.68,
             r"$\alpha_{12}=\alpha_{13}=-\lambda_1=-\frac{\lambda_2}{3}$"
             "\n"
             r"$\Rightarrow \lambda_2 = 3\lambda_1$")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.88, 0.75)
    ax.axis("off")
    ax.set_title(r"(c) Isosceles, $M=3$", fontsize=11, pad=4)

    # ── (d) Equilateral θ=0 ───────────────────────────────────────────────
    ax = fig.add_axes([right_x, bot_bottom, col_w, bot_h])
    a_eq = 0.5
    h_eq = 1 / np.sqrt(3)
    y1 = (0, h_eq)
    y2 = (-a_eq, -h_eq / 2)
    y3 = (a_eq, -h_eq / 2)
    draw_edge(ax, y1, y2, r"$\alpha_{12}=-\frac{\lambda_2}{3}$",
              label_offset_scale=0.16, label_fs=6.5)
    draw_edge(ax, y1, y3, r"$\alpha_{13}=-\frac{\lambda_2}{3}$",
              label_offset_scale=0.16, label_fs=6.5)
    draw_edge(ax, y2, y3, r"$\alpha_{23}=\frac{\lambda_2/3 - \lambda_1}{2}$",
              label_offset_scale=0.10, label_fs=6.5)
    draw_peak(ax, y1, r"$y_1$")
    draw_peak(ax, y2, r"$y_2$", offset=(-0.12, 0.08))
    draw_peak(ax, y3, r"$y_3$")
    ev_origin = (-0.82, -0.48)
    draw_eigvec(ax, ev_origin, (0.33, 0), r"$\mathbf{r}_1$", "#c0392b")
    draw_eigvec(ax, ev_origin, (0, 0.33), r"$\mathbf{r}_2$", "#2980b9")
    info_box(ax, 0, -0.68,
             r"$a = \frac{1}{2},\; h = \frac{1}{\sqrt{3}}$"
             "\n"
             r"$\alpha_{23} \geq 0 \Leftrightarrow |\lambda_2|/|\lambda_1| \geq 3$")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.88, 0.75)
    ax.axis("off")
    ax.set_title(r"(d) Equilateral, $\theta = 0$", fontsize=11, pad=4)

    save_figure(fig, OUTPUT_DIR, "figure_2_peak_configurations")

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE 2B: Panel (e) — Perturbation δy⁰ at R=0
    # ══════════════════════════════════════════════════════════════════════════
    fig_e, ax = plt.subplots(figsize=(5.5, 4.8))

    mu_val = 4.0
    theta_rot = np.radians(19.3)
    ct, st = np.cos(theta_rot), np.sin(theta_rot)
    R_mat = np.array([[ct, -st], [st, ct]])
    h_eq = 1 / np.sqrt(3)
    a_eq = 0.5
    y1_base = np.array([0, h_eq])
    y2_base = np.array([-a_eq, -h_eq / 2])
    y3_base = np.array([a_eq, -h_eq / 2])
    y1_arr = R_mat @ y1_base
    y2_arr = R_mat @ y2_base
    y3_arr = R_mat @ y3_base

    a12_v, a13_v, a23_v = compute_alphas(mu_val, theta_rot)
    beta = 0.5 * np.log(a12_v / a23_v)
    gamma_v = 2 * beta / (3 * np.sqrt(3))
    r1, r2 = np.array([1, 0]), np.array([0, 1])
    dy1 = -gamma_v * r2
    dy2 = -beta * r1 + (gamma_v / 2) * r2
    dy3 = beta * r1 + (gamma_v / 2) * r2
    sc = 0.25
    dy1s, dy2s, dy3s = dy1 * sc, dy2 * sc, dy3 * sc

    # Dashed original triangle
    for p1, p2 in [(y1_arr, y2_arr), (y1_arr, y3_arr), (y2_arr, y3_arr)]:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "--", color="#bdc3c7", lw=1.0, zorder=1)
    # Solid perturbed triangle
    py1, py2, py3 = y1_arr + dy1s, y2_arr + dy2s, y3_arr + dy3s
    for p1, p2 in [(py1, py2), (py1, py3), (py2, py3)]:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color="#2c3e50", lw=1.5, zorder=2)
    # Unperturbed markers
    for pos in [y1_arr, y2_arr, y3_arr]:
        ax.plot(*pos, "o", color="#bdc3c7", ms=8, zorder=3)
    # Perturbed markers
    for pos, col in [(py1, "#e74c3c"), (py2, "#27ae60"), (py3, "#2980b9")]:
        ax.plot(*pos, "s", color=col, ms=9, zorder=4)
    # Perturbation arrows
    for orig, pert, col in [
        (y1_arr, py1, "#e74c3c"),
        (y2_arr, py2, "#27ae60"),
        (y3_arr, py3, "#2980b9"),
    ]:
        ax.annotate(
            "", xy=pert, xytext=orig,
            arrowprops=dict(arrowstyle="->", color=col, lw=1.5),
        )

    # Eigenvectors — bottom-center, well clear of geometry
    ev_origin = (-0.15, -0.70)
    draw_eigvec(ax, ev_origin, (0.30, 0), r"$\mathbf{r}_1$", "#c0392b", fs=9)
    draw_eigvec(ax, ev_origin, (0, 0.30), r"$\mathbf{r}_2$", "#2980b9", fs=9)

    # Formula box
    ax.text(
        0.15, -0.92,
        r"$\delta y_1^0 = -\gamma\mathbf{r}_2$" ", "
        r"$\delta y_{2,3}^0 = \mp\beta\mathbf{r}_1 + \frac{\gamma}{2}\mathbf{r}_2$"
        "\n"
        r"$\beta = \frac{1}{2}\ln(\alpha_{12}/\alpha_{23})$" ", "
        r"$\gamma = \frac{2\beta}{3\sqrt{3}}$",
        fontsize=8.5, ha="center", va="top", color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#e8f8f5", edgecolor="#a3e4d7"),
    )

    # Legend
    ax.plot([], [], "o", color="#bdc3c7", ms=7, label=r"$y_i$ (unperturbed)")
    ax.plot([], [], "s", color="#555", ms=7, label=r"$y_i + \delta y_i^0$ (Prop. 4.1)")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.9, borderpad=0.3,
              handletextpad=0.4)

    ax.set_xlim(-0.90, 1.05)
    ax.set_ylim(-1.18, 0.72)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        r"(e) Perturbation $\delta y^0$ at $R=0$, $\mu=4$, $\theta \approx 19^\circ$",
        fontsize=11, pad=5,
    )

    save_figure(fig_e, OUTPUT_DIR, "figure_2b_perturbation")
    print("Done!")


def create_figure_3():
    """Figure 3 — α vs θ for representative μ values (2×2 subplots)."""
    print(
        "Creating Figure 3 - α vs θ for representative μ values...", end=" ", flush=True
    )

    theta = np.linspace(0, np.pi / 2, 500)
    mu_vals = [0.3, 1.0, 3.0, 6.0]
    mu_colors = ["#8e44ad", "#2980b9", "#27ae60", "#e67e22"]

    # Consistent coefficient styling across all panels
    coeff_styles = {
        "a12": ("-",  2.0, "#c0392b", r"$\alpha_{12}$"),
        "a13": ("--", 1.8, "#2980b9", r"$\alpha_{13}$"),
        "a23": (":",  2.2, "#27ae60", r"$\alpha_{23}$"),
    }

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.28, wspace=0.18)

    for ax, mu, panel_col in zip(axes.flat, mu_vals, mu_colors):
        a12, a13, a23 = compute_alphas(mu, theta)

        # Forbidden region shading
        ax.axhspan(-4, 0, color="#c0392b", alpha=0.05, zorder=0)
        ax.axhline(0, color="k", lw=0.6, zorder=1)

        for data, key in [(a12, "a12"), (a13, "a13"), (a23, "a23")]:
            ls, lw, col, label = coeff_styles[key]
            ax.plot(
                np.degrees(theta),
                data,
                ls=ls,
                lw=lw,
                color=col,
                label=label,
                zorder=3,
            )

        # Check if any coefficient goes negative
        all_positive = np.all(a12 >= -1e-12) and np.all(a13 >= -1e-12) and np.all(a23 >= -1e-12)
        if not all_positive:
            # Find where each goes negative and shade lightly
            for data, key in [(a12, "a12"), (a13, "a13"), (a23, "a23")]:
                neg_mask = data < -1e-12
                if np.any(neg_mask):
                    theta_deg = np.degrees(theta)
                    ax.fill_between(
                        theta_deg, -4, 0,
                        where=neg_mask,
                        color=coeff_styles[key][2],
                        alpha=0.06,
                        zorder=0,
                    )

        ax.set_xlim(0, 90)
        ax.set_ylim(-0.8, 2.5)
        ax.set_title(rf"$\mu = {mu}$", fontsize=12, pad=8, color=panel_col, fontweight="bold")
        ax.grid(True, alpha=0.15, lw=0.5)

        # Admissibility annotation
        if all_positive:
            ax.text(
                0.5, 0.02, r"All $\theta$ admissible",
                transform=ax.transAxes, fontsize=8, ha="center",
                color="#27ae60", style="italic",
            )

    # Axis labels on outer edges only
    for ax in axes[1, :]:
        ax.set_xlabel(r"$\theta$ (degrees)", labelpad=6)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\alpha_{ij}(\theta)$", labelpad=8)

    # Single shared legend
    handles = [
        Line2D([0], [0], color=v[2], ls=v[0], lw=v[1], label=v[3])
        for v in coeff_styles.values()
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        fontsize=9,
        framealpha=0.95,
        edgecolor="#ccc",
        bbox_to_anchor=(0.5, -0.01),
    )

    fig.tight_layout(pad=1.2, rect=[0, 0.04, 1, 1])
    save_figure(fig, OUTPUT_DIR, "figure_3_alpha_vs_theta_vs_mu")
    print("Done!")


def create_figure_4():
    """Figure 4 — Multi-Peak Profile Decomposition Illustration."""
    print(
        "Creating Figure 4 - Profile Decomposition Illustration...", end=" ", flush=True
    )

    def sech_profile(x, z, amp=1.0, width=1.0):
        return amp / np.cosh((x - z) / width)

    x = np.linspace(-14, 14, 1200)
    z1, z2, z3 = -5.0, 0.5, 6.0
    u1 = sech_profile(x, z1, 1.0, 1.0)
    u2 = sech_profile(x, z2, 0.95, 1.0)
    u3 = sech_profile(x, z3, 1.0, 1.0)
    u_sum = u1 + u2 + u3
    h_R = 0.035 * np.sin(1.8 * x) * np.exp(-0.015 * x**2)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 6), height_ratios=[3, 1], gridspec_kw={"hspace": 0.32}
    )

    ax1.fill_between(x, 0, u_sum, color="#4a90d9", alpha=0.12)
    ax1.plot(
        x,
        u_sum,
        color="#1a3a6e",
        lw=2.2,
        label=r"$u_E(x) = \sum_{i=1}^{M} u_R(x - z_i) + h_R(x)$",
    )
    ax1.plot(x, u1, "--", color="#e74c3c", lw=1.1, label=r"$u_R(x - z_1)$")
    ax1.plot(x, u2, "--", color="#27ae60", lw=1.1, label=r"$u_R(x - z_2)$")
    ax1.plot(x, u3, "--", color="#f39c12", lw=1.1, label=r"$u_R(x - z_3)$")

    for zi, lab, col in [
        (z1, r"$z_1$", "#e74c3c"),
        (z2, r"$z_2$", "#27ae60"),
        (z3, r"$z_3$", "#f39c12"),
    ]:
        ax1.axvline(zi, color=col, ls=":", lw=0.6, alpha=0.4)
        ax1.text(zi, -0.12, lab, ha="center", fontsize=10, color=col)

    # Minimum inter-peak distance — shifted down to avoid peaks
    ax1.annotate(
        "",
        xy=(z2, 0.42),
        xytext=(z1, 0.42),
        arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9),
    )
    ax1.text(
        (z1 + z2) / 2,
        0.48,
        r"$m_R = \min_{i \neq j} |z_i - z_j|$",
        ha="center",
        fontsize=9,
        color="#555",
        bbox=dict(
            boxstyle="round,pad=0.15", facecolor="white", edgecolor="#ddd", alpha=0.85
        ),
    )

    ax1.text(
        0.02,
        0.95,
        r"Each $u_R$ solves $(-\Delta + R^2 V(x_0) + 1)u + \sigma|u|^{2p}u = 0$",
        transform=ax1.transAxes,
        fontsize=8,
        va="top",
        color="#555",
        style="italic",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#fef9e7", edgecolor="#ddd"),
    )

    ax1.set_ylabel(r"$u_E(x)$", labelpad=8)
    ax1.set_xlim(-14, 14)
    ax1.set_ylim(-0.18, 1.25)
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9, borderpad=0.7)
    ax1.set_title(
        r"Multi-Peak Decomposition: $u_E(x) \approx \sum_{i=1}^{M} u_R(x - z_i(R))$",
        fontsize=12,
        pad=10,
    )

    ax2.plot(x, h_R, color="#8e44ad", lw=1.2)
    ax2.axhline(0, color="k", lw=0.3)
    ax2.fill_between(x, 0, h_R, color="#8e44ad", alpha=0.12)
    ax2.set_xlabel(r"$x$", labelpad=6)
    ax2.set_ylabel(r"$h_R(x)$", labelpad=8)
    ax2.set_xlim(-14, 14)
    ax2.set_title(
        r"Remainder $h_R$: $\|h_R\|_{H^1} \to 0$ as $E \to \infty$ (solved via IFT, Lemma 2)",
        fontsize=10,
        pad=8,
    )

    save_figure(fig, OUTPUT_DIR, "figure_4_profile_decomposition")
    print("Done!")


def create_figure_5():
    """Figure 5 — Concentration Points vs Peak Positions (Claim 2 illustration)."""
    print(
        "Creating Figure 5 - Concentration Points vs Peak Positions...",
        end=" ",
        flush=True,
    )

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.2], hspace=0.40, wspace=0.38)

    ax_profiles = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    ax_conv = fig.add_subplot(gs[:, 2])

    x = np.linspace(-10, 10, 1000)
    x0 = 0.0

    def sech(x, z, w):
        return 1.0 / np.cosh((x - z) / w)

    configs = [
        {
            "E": 1,
            "R": 1.0,
            "peaks": [(0, 1.0, 1.5)],
            "title": r"$E = 1$, $R = 1.0$",
            "label": "Single peak at $x_0$",
        },
        {
            "E": 4,
            "R": 0.5,
            "peaks": [(-1.5, 0.9, 0.8), (1.5, 0.9, 0.8)],
            "title": r"$E = 4$, $R = 0.5$",
            "label": "Peaks splitting from $x_0$",
        },
        {
            "E": 16,
            "R": 0.25,
            "peaks": [(-3.0, 1.0, 0.5), (3.0, 1.0, 0.5)],
            "title": r"$E = 16$, $R = 0.25$",
            "label": "Well-separated peaks",
        },
        {
            "E": 64,
            "R": 0.125,
            "peaks": [(-4.0, 1.0, 0.35), (0.0, 0.8, 0.35), (4.0, 1.0, 0.35)],
            "title": r"$E = 64$, $R = 0.125$",
            "label": "$M=3$ peaks",
        },
    ]

    for idx, (ax, cfg) in enumerate(zip(ax_profiles, configs)):
        psi = np.zeros_like(x)
        for z, _, w in cfg["peaks"]:
            psi += sech(x, z, w)

        ax.fill_between(x, 0, psi, color="#4a90d9", alpha=0.15)
        ax.plot(x, psi, color="#1a3a6e", lw=2.0)
        ax.axvline(x0, color="#27ae60", ls="--", lw=1.2, alpha=0.7)
        ax.text(x0 + 0.3, max(psi) * 0.90, r"$x_0$", fontsize=10, color="#27ae60")

        for i, (z, _, w) in enumerate(cfg["peaks"]):
            ax.plot(z, sech(z, z, w), "v", color="#e74c3c", ms=8, zorder=5)
            ax.text(z, -0.15, rf"$z_{i + 1}$", ha="center", fontsize=9, color="#e74c3c")
            Rz = cfg["R"] * z
            ax.plot(Rz, -0.05, "s", color="#8e44ad", ms=5, zorder=5, clip_on=False)

        ax.axhline(0, color="k", lw=0.3)
        ax.set_xlim(-8, 8)
        ax.set_ylim(-0.22, 1.3)
        ax.set_title(cfg["title"], fontsize=10, pad=6)
        if idx >= 2:
            ax.set_xlabel(r"$x$", labelpad=4)
        ax.set_ylabel(r"$|\psi_E(x)|$", labelpad=6)

        ax.text(
            0.97,
            0.97,
            cfg["label"],
            transform=ax.transAxes,
            fontsize=7,
            ha="right",
            va="top",
            color="#555",
            style="italic",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="#f8f8f8", edgecolor="#ddd"),
        )

    ax_profiles[0].plot(
        [], [], "v", color="#e74c3c", ms=7, label=r"Peak centers $z_i(E)$"
    )
    ax_profiles[0].plot([], [], "s", color="#8e44ad", ms=5, label=r"Rescaled $Rz_i$")
    ax_profiles[0].plot(
        [], [], "--", color="#27ae60", lw=1.2, label=r"Critical point $x_0$"
    )
    ax_profiles[0].legend(fontsize=7, loc="upper left", framealpha=0.9, borderpad=0.6)

    # Convergence plot
    E_vals = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    R_vals = E_vals ** (-0.5)
    np.random.seed(42)
    c1, c2 = 3.5, -3.5
    noise_scale = 0.02
    Rz_1 = x0 + c1 * R_vals * (1 + noise_scale * np.random.randn(len(R_vals)))
    Rz_2 = x0 + c2 * R_vals * (1 + noise_scale * np.random.randn(len(R_vals)))

    ax_conv.semilogx(
        E_vals, Rz_1, "o-", color="#e74c3c", lw=1.5, ms=6, label=r"$Rz_1(E)$"
    )
    ax_conv.semilogx(
        E_vals, Rz_2, "s-", color="#2980b9", lw=1.5, ms=6, label=r"$Rz_2(E)$"
    )
    ax_conv.axhline(
        x0, color="#27ae60", ls="--", lw=1.5, label=r"$x_0$ (critical point)"
    )
    ax_conv.fill_between(E_vals, x0 - 0.05, x0 + 0.05, color="#27ae60", alpha=0.1)

    ax_conv.set_xlabel(r"Energy $E$", fontsize=12, labelpad=6)
    ax_conv.set_ylabel(r"Rescaled position $Rz_i(E)$", fontsize=11, labelpad=8)
    ax_conv.set_title(
        r"$Rz_i(E) \to x_0$ as $E \to \infty$" "\n(Proposition 2.3)", fontsize=12, pad=8
    )
    ax_conv.legend(fontsize=9, loc="upper right", borderpad=0.6)
    ax_conv.set_ylim(-2.5, 2.5)
    ax_conv.grid(True, alpha=0.2, which="both")
    ax_conv.text(
        100,
        0.5,
        r"$|Rz_i - x_0| = O(E^{-1/2})$",
        fontsize=10,
        color="#555",
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )

    save_figure(fig, OUTPUT_DIR, "figure_5_concentration_vs_peaks")
    print("Done!")


def create_figure_6():
    """Figure 6 — Interaction Functional Decay."""
    print("Creating Figure 6 - Interaction Functional Decay...", end=" ", flush=True)

    d = np.linspace(0.5, 12, 500)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    R_vals = [0.01, 0.1, 0.3, 0.5, 0.8]
    V0 = 1.0
    colors = cm.get_cmap("plasma")(np.linspace(0.15, 0.85, len(R_vals)))

    for R, col in zip(R_vals, colors):
        kappa = np.sqrt(1 + R**2 * V0)
        Q = np.exp(-d * kappa) * d ** (-0.5)
        ax1.semilogy(d, Q, color=col, lw=1.8, label=rf"$R = {R}$")

    ax1.set_xlabel(r"Inter-peak distance $d = |z_k - z_i|$", labelpad=6)
    ax1.set_ylabel(r"$Q_R(d)$ (log scale)", labelpad=8)
    ax1.legend(fontsize=9, title=r"$R = E^{-1/2}$", title_fontsize=9, borderpad=0.6)
    ax1.set_xlim(0.5, 12)
    ax1.set_ylim(1e-8, 1)
    ax1.set_title(r"Exponential Decay of $Q_R(d)$", fontsize=12, pad=10)

    ax1.text(
        3.5,
        2e-2,
        r"$Q_R(d) \sim C \cdot d^{-(n-1)/2} \cdot e^{-d\sqrt{1 + R^2 V(x_0)}}$",
        fontsize=9,
        color="#333",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax1.axhspan(0, 1e-6, color="#f0f0f0", alpha=0.6, zorder=0)
    ax1.text(
        8.5,
        3e-7,
        "Negligible interaction\n" r"(justifies $N_i$ truncation)",
        fontsize=8,
        color="#7f8c8d",
        ha="center",
        style="italic",
    )

    # Right panel
    d_ratio = np.linspace(1, 3.5, 500)
    for R, col in zip(R_vals, colors):
        kappa = np.sqrt(1 + R**2 * V0)
        m_R = 3.0
        alpha = np.exp(-m_R * kappa * (d_ratio - 1)) * d_ratio ** (-0.5)
        ax2.plot(d_ratio, alpha, color=col, lw=1.8, label=rf"$R = {R}$")

    ax2.axhline(1, color="k", ls=":", lw=0.5)
    ax2.axhline(0, color="k", lw=0.3)
    ax2.set_xlabel(r"Distance ratio $|z_k - z_i| \, / \, m_R$", labelpad=6)
    ax2.set_ylabel(r"$\alpha_{ki} = Q_R(|z_k - z_i|) \, / \, Q_R(m_R)$", labelpad=8)
    ax2.set_xlim(1, 3.5)
    ax2.set_ylim(-0.05, 1.15)
    ax2.set_title(r"Interaction Coefficient $\alpha_{ki}$", fontsize=12, pad=10)

    ax2.annotate(
        "Nearest neighbors\n" r"$|z_k - z_i| = m_R$" "\n" r"$\alpha_{ki} = 1$",
        xy=(1.02, 0.98),
        xytext=(1.4, 0.62),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="#555", lw=0.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#d5f5e3", edgecolor="#82e0aa"),
    )
    ax2.annotate(
        "Next-nearest neighbors\n" r"$\alpha_{ki} \ll 1$" "\n(dropped from sum)",
        xy=(2.0, 0.03),
        xytext=(2.3, 0.45),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="#555", lw=0.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#fadbd8", edgecolor="#f1948a"),
    )
    ax2.text(
        0.02,
        0.02,
        r"$0 \leq \alpha_{ki} \leq 1$ by construction"
        "\n"
        r"$\alpha_{ki} = 1$ iff $k \in N_i$",
        transform=ax2.transAxes,
        fontsize=8,
        va="bottom",
        color="#555",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax2.legend(fontsize=8, loc="center right", borderpad=0.6)

    fig.tight_layout(pad=1.5, w_pad=3.0)
    save_figure(fig, OUTPUT_DIR, "figure_6_interaction_decay")
    print("Done!")


def create_figure_7():
    """Figure 7 — Analytical Pipeline: PDE → Multi-Peak Ground States"""
    print("Creating Figure 7 – Methodology Flowchart...", end=" ", flush=True)

    plt.rcParams.update({"font.size": 7.5})

    def box(ax, x, y, w, h, text, fc, ec, fs=7, bold=False, lw=1.2, ls="-"):
        b = FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.08",
            facecolor=fc,
            edgecolor=ec,
            lw=lw,
            linestyle=ls,
        )
        ax.add_patch(b)
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=fs,
            fontweight="bold" if bold else "normal",
            multialignment="center",
            color="#1a1a1a",
            linespacing=1.25,
        )

    def arr(ax, x1, y1, x2, y2, color=C_ARROW, lw=1.0, ls="-"):
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="->", color=color, lw=lw, linestyle=ls, shrinkA=1, shrinkB=1
            ),
        )

    def col_header(ax, x, y, w, text, color):
        rect = FancyBboxPatch(
            (x - w / 2, y - 0.22),
            w,
            0.44,
            boxstyle="round,pad=0.06",
            facecolor=color,
            edgecolor=color,
            lw=1.5,
            alpha=0.90,
        )
        ax.add_patch(rect)
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=7.5,
            fontweight="bold",
            color="white",
        )

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.set_xlim(-0.2, 12.2)
    ax.axis("off")
    ax.set_aspect("equal")

    Lx = 2.75
    Rx = 9.25
    Bw = 4.6
    sp = 1.12  # uniform step — both columns use this exclusively
    bh = 0.62  # standard box height
    bhL = 0.74  # tall box height  (3-line content)
    tag_h = 0.32  # narrow tag/badge height
    pad = 0.18  # border padding

    Y_top = 17.0

    col_header(ax, Lx, Y_top, Bw + 0.2, "Kirr\u2013Natarajan (2018)", C_BLUE_FG)
    col_header(ax, Rx, Y_top, Bw + 0.2, "Kirr\u2013Knox (2026)", C_HDR_R)

    # ─── helper: step down and draw arrow from previous box bottom ──────────
    # prev_h is the height of the box we're leaving; next_h is the box we arrive at
    def step(ax, cx, Y_prev, prev_h, next_h):
        Y_new = Y_prev - sp
        arr(ax, cx, Y_prev - prev_h / 2, cx, Y_new + next_h / 2)
        return Y_new

    # ═══════════════ LEFT COLUMN ════════════════════════════════════════════
    # 8 boxes on a strict sp grid:
    #   1 PDE (bh)  2 Rescaling (bh)  3 Conc-Comp (bh)  4 Lyap-Schm (bh)
    #   5 Reduced (bhL)  6 Thm 4.3 (bhL)  7 Resolved tag (tag_h)  8 Open (bhL)

    Y = Y_top - sp * 0.71  # first box center  (same as original Y_top - 0.80)

    box(
        ax,
        Lx,
        Y,
        Bw,
        bh,
        r"$(-\Delta + V + E)\psi + \sigma|\psi|^{2p}\psi = 0$"
        "\nFull PDE on "
        r"$\mathbb{R}^n$",
        C_BLUE_BG,
        C_BLUE_FG,
        fs=6.8,
        bold=True,
    )

    Y = step(ax, Lx, Y, bh, bh)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bh,
        r"Rescaling:  $u_E = E^{-1/(2p)}\psi_E(E^{-1/2}\cdot)$"
        "\n"
        r"$\|u_E\|_{H^1}$ bounded;  $R = E^{-1/2}\to 0$",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.3,
    )

    Y = step(ax, Lx, Y, bh, bh)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bh,
        "Conc.-Compactness (Lions)"
        "\n"
        r"$u_E \approx \sum_{i=1}^{M} u_R(\cdot - z_i) + h_R$",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.3,
    )

    Y = step(ax, Lx, Y, bh, bh)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bh,
        "Lyapunov\u2013Schmidt (Lemma 2)"
        "\n"
        r"$P_\perp F = 0 \Rightarrow h_R = h_R(z_1,\ldots,z_M,R)$",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.3,
    )

    Y = step(ax, Lx, Y, bh, bhL)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bhL,
        "Reduced System"
        "\n"
        r"$H_V y_i = \sum_{k \in N_i} \alpha_{ki}\,\chi(y_k - y_i)$"
        "\n"
        r"peaks $\in$ unstable manifold of $H_V$",
        C_GREEN_BG,
        C_GREEN_FG,
        fs=6.3,
        bold=True,
    )

    Y = step(ax, Lx, Y, bhL, bhL)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bhL,
        "Thm 4.3: one profile per critical point"
        "\n"
        r"$\|U_E - \sum \mu_i u_\infty(\cdot - x_i\sqrt{E})\|_{H^2} < \varepsilon$"
        "\nExistence + uniqueness + stability",
        C_GREEN_BG,
        C_GREEN_FG,
        fs=6.3,
    )

    Y = step(ax, Lx, Y, bhL, tag_h)
    box(
        ax,
        Lx,
        Y,
        Bw * 0.78,
        tag_h,
        "Resolved for distinct critical points",
        C_GREEN_BG,
        "#1a5632",
        fs=6.0,
        bold=True,
    )

    Y = step(ax, Lx, Y, tag_h, bhL)
    box(
        ax,
        Lx,
        Y,
        Bw,
        bhL,
        "Open Problem (K-N Remark 4.1):"
        "\n"
        r"Multiple profiles $\to$ same critical point"
        "\n"
        r"along eigenvectors with $\lambda_j < 0$",
        C_RED_BG,
        C_RED_FG,
        fs=6.3,
        bold=True,
        lw=1.6,
    )

    Y_open = Y
    yb_L = Y_open - bhL / 2 - 0.10
    yt_col = Y_top - 0.48

    left_rect = FancyBboxPatch(
        (Lx - Bw / 2 - pad, yb_L),
        Bw + 2 * pad,
        yt_col - yb_L,
        boxstyle="round,pad=0.10",
        facecolor="none",
        edgecolor=C_BLUE_FG,
        lw=0.7,
        alpha=0.22,
    )
    ax.add_patch(left_rect)

    # ═══════════════ BRIDGE ═════════════════════════════════════════════════
    Y_H = Y_top - sp * 0.71  # top of right column (same as left)
    mid_x = (Lx + Bw / 2 + Rx - Bw / 2) / 2
    mid_y = (Y_open + Y_H) / 2
    arr(
        ax,
        Lx + Bw / 2 + 0.08,
        Y_open + 0.05,
        Rx - Bw / 2 - 0.08,
        Y_H - 0.05,
        color=C_RED_FG,
        lw=1.5,
    )
    # ax.text(mid_x + 0.1, mid_y + 0.15, 'this work', ha='center', va='bottom', fontsize=5.5, fontstyle='italic', color=C_RED_FG)

    # ═══════════════ RIGHT COLUMN ═══════════════════════════════════════════
    # 8 boxes on the same sp grid:
    #   1 Classify (bh)  2 Admissibility (bh)  3 Perturbation (bh)
    #   4 Explicit Sol (bhL)  5 Kernel Obstruction (bhL)  6 Decompose (bhL)
    #   7 Gauge-Fixing (bh)  8 IFT+Lift (bhL)

    Y = Y_H

    box(
        ax,
        Rx,
        Y,
        Bw,
        bh,
        r"Classify Limiting Configs ($R = 0$)"
        "\nCollinear | Isosceles | Equilateral | Rot. "
        r"$\theta$",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.3,
    )

    Y = step(ax, Rx, Y, bh, bh)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bh,
        r"Admissibility:  $\alpha_{ki} \geq 0$"
        "\n"
        r"$\frac{1}{3} \leq \mu \leq 3$: all $\theta$;  restricted otherwise",
        C_GREEN_BG,
        C_GREEN_FG,
        fs=6.3,
    )

    Y = step(ax, Rx, Y, bh, bh)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bh,
        r"Perturbation System:  $F_i(\delta y, R) = 0$"
        "\n"
        r"Extended to $R = 0$ by continuity",
        C_BLUE_BG,
        C_BLUE_FG,
        fs=6.3,
    )

    Y = step(ax, Rx, Y, bh, bhL)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bhL,
        r"Explicit Solution at $R=0$ (Claim 3)"
        "\n"
        r"$\delta y_1^0=-\gamma\mathbf{r}_2$, "
        r"$\delta y_{2,3}^0=\mp\beta\mathbf{r}_1+\frac{\gamma}{2}\mathbf{r}_2$",
        C_BLUE_BG,
        C_BLUE_FG,
        fs=6.3,
    )

    Y = step(ax, Rx, Y, bhL, bhL)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bhL,
        "Kernel Obstruction (Claim 4)"
        "\n"
        r"$\Omega = \ker(D_{\delta y}F_i) \supseteq \{(v,v,v):v \in \mathbb{R}^2\}$"
        "\nIFT blocked by translational invariance",
        C_RED_BG,
        C_RED_FG,
        fs=6.3,
        bold=True,
        lw=1.6,
    )

    Y = step(ax, Rx, Y, bhL, bhL)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bhL,
        r"Decompose:  $\mathbb{R}^{2\times 3}=\Omega\oplus\Omega^\perp$"
        "\n"
        r"$P_\perp F$: IFT on $\Omega^\perp$  |  "
        r"$P_w F$: via equivariance",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.0,
    )

    Y = step(ax, Rx, Y, bhL, bh)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bh,
        r"Gauge-Fixing:  $\delta y_1 \perp \mathbf{r}_1$  ;  "
        r"$(\delta y_3-\delta y_2) \parallel \mathbf{r}_1$",
        C_GRAY_BG,
        C_GRAY_FG,
        fs=6.3,
    )

    Y = step(ax, Rx, Y, bh, bhL)
    box(
        ax,
        Rx,
        Y,
        Bw,
        bhL,
        "IFT + Lift via Lemma 2"
        "\n"
        r"$\Rightarrow$ exact multi-peak $\psi_E$ for $E \gg 1$"
        "\n(contingent on completing open steps)",
        C_GREEN_BG,
        "#1a5632",
        fs=6.3,
        bold=True,
        lw=1.3,
        ls="--",
    )
    Y_Q = Y

    yb_R = Y_Q - bhL / 2 - 0.10
    right_rect = FancyBboxPatch(
        (Rx - Bw / 2 - pad, yb_R),
        Bw + 2 * pad,
        yt_col - yb_R,
        boxstyle="round,pad=0.10",
        facecolor="none",
        edgecolor=C_HDR_R,
        lw=0.7,
        alpha=0.22,
    )
    ax.add_patch(right_rect)

    # Open-steps callout (below right border, same as original)
    Y_os = yb_R - 0.42
    box(
        ax,
        Rx,
        Y_os,
        Bw,
        0.55,
        "Remaining:  (1) Verify equivariance of "
        r"$F_i$"
        "\n(2) Invertibility on gauge-fixed subspace"
        "\n(3) Extend to continuous "
        r"$\theta$",
        C_OPEN_BG,
        C_OPEN_FG,
        fs=5.8,
        lw=1.0,
        ls="--",
    )
    arr(ax, Rx, Y_os + 0.28, Rx, yb_R, color=C_OPEN_FG, lw=0.7, ls="--")

    # ═══════════════ LEGEND ═════════════════════════════════════════════════
    # Anchor to fixed position just below the lower of the two column bottoms
    # ═══════════════ LEGEND ═════════════════════════════════════════════════
    leg_y = min(yb_L, Y_os - 0.55 / 2) - 0.52
    dx = 2.20  # was 3.05 — tighter spacing between items
    total_w = 3 * dx + 0.60  # 3 gaps + width of last item
    leg_x = (12.0 - total_w) / 2  # center the group across the 12-unit axis
    for i, (fc, ec, label) in enumerate(
        [
            (C_BLUE_BG, C_BLUE_FG, "Systems"),
            (C_GRAY_BG, C_GRAY_FG, "Analysis"),
            (C_GREEN_BG, C_GREEN_FG, "Results"),
            (C_RED_BG, C_RED_FG, "Open"),
        ]
    ):
        xx = leg_x + i * dx
        swatch = FancyBboxPatch(
            (xx, leg_y - 0.09),
            0.26,
            0.18,
            boxstyle="round,pad=0.03",
            facecolor=fc,
            edgecolor=ec,
            lw=0.7,
        )
        ax.add_patch(swatch)
        ax.text(xx + 0.34, leg_y, label, fontsize=5.5, va="center", color="#333")

    ax.set_ylim(leg_y - 0.40, 17.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    save_figure(fig, OUTPUT_DIR, "figure_7_methodology_flowchart", pad_inches=0.08)
    print("Done!")


if __name__ == "__main__":
    figures = [
        create_figure_1,
        create_figure_2,
        create_figure_3,
        create_figure_4,
        create_figure_5,
        create_figure_6,
        create_figure_7,
    ]
    for create in figures:
        create()
