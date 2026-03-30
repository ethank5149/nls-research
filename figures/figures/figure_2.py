import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from config import OUTPUT_DIR

warnings.filterwarnings("ignore")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",  # type: ignore[dict-item]
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.5,
    }
)


def compute_alphas(mu, theta):
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


def box(ax, x, y, w, h, text, fc, ec, fs=9, bold=False):
    b = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.15",
        facecolor=fc,
        edgecolor=ec,
        lw=1.5,
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


def arr(ax, x1, y1, x2, y2, color="#2c3e50", lw=1.5):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", color=color, lw=lw),
    )


def create():
    print("Creating Figure 2 - Extended Peak Configurations...", end=" ", flush=True)
    fig = plt.figure(figsize=(16, 9))

    # Grid: 2 rows. Top row: 3 panels. Bottom row: 2 panels centered.
    ax_a = fig.add_subplot(231)
    ax_b = fig.add_subplot(232)
    ax_c = fig.add_subplot(233)
    ax_d = fig.add_subplot(234)
    ax_e = fig.add_subplot(235)
    # Leave 236 empty for balance

    def draw_peak(ax, pos, label, color="#2c3e50", ms=9):
        ax.plot(*pos, "o", color=color, ms=ms, zorder=5)
        ax.annotate(
            label,
            xy=pos,
            xytext=(pos[0] + 0.06, pos[1] + 0.1),
            fontsize=10,
            color=color,
            ha="center",
        )  # type: ignore[arg-type]

    def draw_edge(ax, p1, p2, alpha_val=None, color="#7f8c8d", lw=1.3):
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color=color, lw=lw, zorder=2)
        if alpha_val is not None:
            mx, my = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
            # Offset perpendicular to edge
            dx, dy = p2[0] - p1[0], p2[1] - p1[1]
            length = np.sqrt(dx**2 + dy**2)
            nx, ny = -dy / length * 0.1, dx / length * 0.1
            ax.text(
                mx + nx,
                my + ny,
                alpha_val,
                fontsize=7,
                ha="center",
                va="center",
                color="#2980b9",
                fontweight="bold",
                bbox=dict(
                    boxstyle="round,pad=0.15",
                    facecolor="white",
                    edgecolor="#ddd",
                    alpha=0.9,
                ),
            )

    def draw_eigvec(ax, origin, direction, label, color="#c0392b"):
        ax.annotate(
            "",
            xy=(origin[0] + direction[0], origin[1] + direction[1]),
            xytext=origin,
            arrowprops=dict(arrowstyle="->", color=color, lw=1.3),
        )  # type: ignore[arg-type]
        ax.text(
            origin[0] + direction[0] * 1.12,
            origin[1] + direction[1] * 1.12,
            label,
            fontsize=9,
            color=color,
            ha="center",
            va="center",
        )

    # ── (a) Collinear M=2 ──
    ax = ax_a
    peaks = [(-0.5, 0), (0.5, 0)]
    draw_edge(ax, peaks[0], peaks[1], r"$\alpha_{12}$")
    draw_peak(ax, peaks[0], r"$y_{-1}$")
    draw_peak(ax, peaks[1], r"$y_1$")
    draw_eigvec(ax, (-0.2, -0.45), (0.5, 0), r"$\mathbf{r}_1$")
    # Distance annotation
    ax.annotate(
        "",
        xy=(0.5, -0.22),
        xytext=(-0.5, -0.22),
        arrowprops=dict(arrowstyle="<->", color="#555", lw=0.8),
    )
    ax.text(0, -0.3, r"$|y_1 - y_{-1}| = 1$", fontsize=8, ha="center", color="#555")
    # Constraint box
    ax.text(
        0,
        -0.55,
        r"$H_V y_{\pm 1} = \alpha_{12}\mathbf{r}_1$"
        "\n"
        r"$\Rightarrow \alpha_{12} = -\lambda_1$",
        fontsize=8,
        ha="center",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.8, 0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(a) Collinear, $M=2$", fontsize=11, pad=8)

    # ── (b) Collinear M=5 ──
    ax = ax_b
    peaks_5 = [(-1, 0), (-0.5, 0), (0, 0), (0.5, 0), (1, 0)]
    labels_5 = [r"$y_{-2}$", r"$y_{-1}$", r"$y_0$", r"$y_1$", r"$y_2$"]
    for i in range(len(peaks_5) - 1):
        draw_edge(ax, peaks_5[i], peaks_5[i + 1])
    for p, l in zip(peaks_5, labels_5):
        draw_peak(ax, p, l)
    draw_eigvec(ax, (-0.3, -0.4), (0.5, 0), r"$\mathbf{r}_1$")
    ax.text(
        0,
        -0.55,
        r"$H_V(i\mathbf{r}_1) = (\alpha_{i,i+1} + \alpha_{i,i-1})\mathbf{r}_1$"
        "\n"
        r"Nearest-neighbor coupling only",
        fontsize=7.5,
        ha="center",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.8, 0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(b) Collinear, $M=5$ (odd)", fontsize=11, pad=8)

    # ── (c) Isosceles triangle ──
    ax = ax_c
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
        draw_edge(ax, p1, p2, albl, color=c_edge)
    draw_peak(ax, y1, r"$y_1$")
    draw_peak(ax, y2, r"$y_2$")
    draw_peak(ax, y3, r"$y_3$")
    draw_eigvec(ax, (-0.15, -0.75), (0.45, 0), r"$\mathbf{r}_1$", "#c0392b")
    draw_eigvec(ax, (-0.75, -0.15), (0, 0.45), r"$\mathbf{r}_2$", "#2980b9")
    ax.text(
        0,
        -0.85,
        r"$\alpha_{12}=\alpha_{13}=-\lambda_1=-\frac{\lambda_2}{3}$"
        "\n"
        r"$\Rightarrow \lambda_2 = 3\lambda_1$",
        fontsize=8,
        ha="center",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    # Geometric constraints
    ax.text(
        0.6,
        0.45,
        r"$a > \frac{1}{2}$" "\n" r"$h < \frac{1}{\sqrt{3}}$",
        fontsize=8,
        color="#8e44ad",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#f5eef8", edgecolor="#d2b4de"),
    )
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.05, 0.85)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(c) Isosceles, $M=3$", fontsize=11, pad=8)

    # ── (d) Equilateral θ=0 with α values ──
    ax = ax_d
    a_eq = 0.5
    h_eq = 1 / np.sqrt(3)
    y1 = (0, h_eq)
    y2 = (-a_eq, -h_eq / 2)
    y3 = (a_eq, -h_eq / 2)
    draw_edge(ax, y1, y2, r"$\alpha_{12}=-\frac{\lambda_2}{3}$")
    draw_edge(ax, y1, y3, r"$\alpha_{13}=-\frac{\lambda_2}{3}$")
    draw_edge(ax, y2, y3, r"$\alpha_{23}=\frac{\lambda_2/3 - \lambda_1}{2}$")
    draw_peak(ax, y1, r"$y_1$")
    draw_peak(ax, y2, r"$y_2$")
    draw_peak(ax, y3, r"$y_3$")
    draw_eigvec(ax, (-0.15, -0.75), (0.45, 0), r"$\mathbf{r}_1$", "#c0392b")
    draw_eigvec(ax, (-0.75, -0.15), (0, 0.45), r"$\mathbf{r}_2$", "#2980b9")
    ax.text(
        0,
        -0.88,
        r"$a = \frac{1}{2},\; h = \frac{1}{\sqrt{3}}$"
        "\n"
        r"$\alpha_{23} \geq 0 \Leftrightarrow |\lambda_2|/|\lambda_1| \geq 3$",
        fontsize=8,
        ha="center",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 0.85)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(d) Equilateral, $\theta = 0$", fontsize=11, pad=8)

    # ── (e) Equilateral with perturbation δy⁰ ──
    ax = ax_e
    y1_arr = np.array([0, h_eq])
    y2_arr = np.array([-a_eq, -h_eq / 2])
    y3_arr = np.array([a_eq, -h_eq / 2])
    mu_val = 4.0
    lam1 = -1.0
    lam2 = -mu_val
    alpha12 = -lam2 / 3
    alpha23_val = (-lam1 + lam2 / 3) / 2
    beta = 0.5 * np.log(alpha12 / alpha23_val)
    gamma = 2 * beta / (3 * np.sqrt(3))
    r1 = np.array([1, 0])
    r2 = np.array([0, 1])
    dy1 = -gamma * r2
    dy2 = -beta * r1 + (gamma / 2) * r2
    dy3 = beta * r1 + (gamma / 2) * r2
    sc = 0.7
    dy1s, dy2s, dy3s = dy1 * sc, dy2 * sc, dy3 * sc

    # Original (dashed)
    for p1, p2 in [(y1_arr, y2_arr), (y1_arr, y3_arr), (y2_arr, y3_arr)]:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "--", color="#bdc3c7", lw=1.0, zorder=1)
    # Perturbed (solid)
    py1, py2, py3 = y1_arr + dy1s, y2_arr + dy2s, y3_arr + dy3s
    for p1, p2 in [(py1, py2), (py1, py3), (py2, py3)]:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], "-", color="#2c3e50", lw=1.5, zorder=2)

    # Original positions (gray)
    for pos in [y1_arr, y2_arr, y3_arr]:
        ax.plot(*pos, "o", color="#bdc3c7", ms=6, zorder=3)
    # Perturbed positions (colored)
    for pos, col in [(py1, "#e74c3c"), (py2, "#27ae60"), (py3, "#2980b9")]:
        ax.plot(*pos, "s", color=col, ms=7, zorder=4)
    # Displacement arrows
    for orig, pert, col in [
        (y1_arr, py1, "#e74c3c"),
        (y2_arr, py2, "#27ae60"),
        (y3_arr, py3, "#2980b9"),
    ]:
        ax.annotate(
            "",
            xy=pert,
            xytext=orig,
            arrowprops=dict(arrowstyle="->", color=col, lw=1.3),
        )  # type: ignore[arg-type]

    draw_eigvec(ax, (-0.15, -0.7), (0.4, 0), r"$\mathbf{r}_1$", "#c0392b")
    draw_eigvec(ax, (-0.7, -0.15), (0, 0.4), r"$\mathbf{r}_2$", "#2980b9")

    ax.text(
        0,
        -0.88,
        r"$\delta y_1^0 = -\gamma\mathbf{r}_2$, "
        "\n"
        r"$\delta y_{2,3}^0 = \mp\beta\mathbf{r}_1 + \frac{\gamma}{2}\mathbf{r}_2$"
        f"\n"
        r"$\beta = \frac{1}{2}\ln(\alpha_{12}/\alpha_{23}),\; "
        "\n"
        r"\gamma = \frac{2\beta}{3\sqrt{3}}$",
        fontsize=7,
        ha="center",
        color="#444",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#e8f8f5", edgecolor="#a3e4d7"),
    )
    # Legend for this panel
    ax.plot([], [], "o", color="#bdc3c7", ms=5, label=r"$y_i$ (unperturbed)")
    ax.plot([], [], "s", color="#555", ms=5, label=r"$y_i + \delta y_i^0$ (Claim 3)")
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)

    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.1, 0.85)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"(e) Perturbation $\delta y^0$ at $R=0$, $\mu=4$", fontsize=11, pad=8)

    # Hide unused subplot
    fig.delaxes(fig.add_subplot(236))

    fig.tight_layout(h_pad=1.5, w_pad=1.0)
    fig.savefig(OUTPUT_DIR / "figure_2_peak_configurations.pdf")
    fig.savefig(OUTPUT_DIR / "figure_2_peak_configurations.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
