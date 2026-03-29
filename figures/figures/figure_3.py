import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


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
    print(
        "Creating Figure 3 - α vs θ for representative μ values...", end=" ", flush=True
    )
    theta = np.linspace(0, np.pi / 2, 500)
    # 4 representative mu values
    mu_vals = [0.3, 1.0, 3.0, 6.0]
    mu_colors = ["#8e44ad", "#2980b9", "#27ae60", "#e67e22"]

    # Linestyles for the three coefficients
    styles = {
        "a12": ("-", 2.0),  # solid, thick
        "a13": ("--", 1.8),  # dashed
        "a23": (":", 2.2),  # dotted, slightly thicker for visibility
    }

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Forbidden region
    ax.axhspan(-4, 0, color="#e74c3c", alpha=0.06, zorder=0)
    ax.axhline(0, color="k", lw=0.6, zorder=1)

    for mu, col in zip(mu_vals, mu_colors):
        a12, a13, a23 = compute_alphas(mu, theta)
        ax.plot(
            np.degrees(theta),
            a12,
            ls=styles["a12"][0],
            lw=styles["a12"][1],
            color=col,
            zorder=3,
        )
        ax.plot(
            np.degrees(theta),
            a13,
            ls=styles["a13"][0],
            lw=styles["a13"][1],
            color=col,
            zorder=3,
        )
        ax.plot(
            np.degrees(theta),
            a23,
            ls=styles["a23"][0],
            lw=styles["a23"][1],
            color=col,
            zorder=3,
        )

    # Forbidden label
    ax.text(
        45,
        -0.25,
        r"$\alpha < 0$: forbidden (focusing nonlinearity requires attractive interaction)",
        fontsize=8,
        ha="center",
        color="#c0392b",
        style="italic",
    )
    ax.set_xlabel(r"$\theta$ (degrees)")
    ax.set_ylabel(r"$\alpha_{ij}(\theta)$")
    ax.set_xlim(0, 90)
    ax.set_ylim(-0.8, 2.5)

    # Two-part legend: linestyle (coefficient type) and color (mu value)
    coeff_handles = [
        Line2D(
            [0],
            [0],
            color="#333",
            ls="-",
            lw=2.0,
            label=r"$\alpha_{12}$  (peak 1 $\leftrightarrow$ 2, apex to base-left)",
        ),
        Line2D(
            [0],
            [0],
            color="#333",
            ls="--",
            lw=1.8,
            label=r"$\alpha_{13}$  (peak 1 $\leftrightarrow$ 3, apex to base-right)",
        ),
        Line2D(
            [0],
            [0],
            color="#333",
            ls=":",
            lw=2.2,
            label=r"$\alpha_{23}$  (peak 2 $\leftrightarrow$ 3, base edge)",
        ),
    ]
    mu_handles = [
        Line2D([0], [0], color=col, ls="-", lw=3, label=rf"$\mu = {mu}$")
        for mu, col in zip(mu_vals, mu_colors)
    ]
    leg1 = ax.legend(
        handles=coeff_handles,
        loc="upper left",
        fontsize=8,
        title="Coefficient (linestyle)",
        title_fontsize=8,
        framealpha=0.95,
        edgecolor="#ccc",
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=mu_handles,
        loc="upper right",
        fontsize=8,
        title=r"$\mu = |\lambda_2|/|\lambda_1|$ (color)",
        title_fontsize=8,
        framealpha=0.95,
        edgecolor="#ccc",
    )
    ax.set_title(
        r"Interaction Coefficients $\alpha_{12}$, $\alpha_{13}$, $\alpha_{23}$ vs Orientation $\theta$",
        fontsize=13,
    )

    fig.tight_layout()
    fig.savefig(f"figures/output/figure_3_alpha_vs_theta_vs_mu.pdf")
    fig.savefig(f"figures/output/figure_3_alpha_vs_theta_vs_mu.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
