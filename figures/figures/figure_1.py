import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import ListedColormap
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

    # Determine binding constraint everywhere
    stack = np.stack([a12, a13, a23], axis=-1)
    binding = np.argmin(stack, axis=-1).astype(float)

    # Mask: inside admissible region show binding constraint, outside show gray
    binding_masked = np.where(admissible, binding, 3.0)  # 3 = inadmissible

    # Aspect ratio to match Mathematica: x-range [0,6], y-range [0,π/2] ≈ [0,1.571]
    # Natural ratio 6 : π/2 ≈ 3.82:1
    fig, ax = plt.subplots(figsize=(6, 6))

    # Custom colormap: 0=alpha12 red, 1=alpha13 orange, 2=alpha23 gray-blue, 3=light gray
    cmap = ListedColormap(["#e74c3c", "#e67e22", "#7f8c8d", "#ececec"])

    ax.pcolormesh(
        MU, TH, binding_masked, cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto"
    )

    # Admissibility boundary (the key contour)
    ax.contour(
        MU, TH, admissible.astype(float), levels=[0.5], colors="#1a1a2e", linewidths=2.0
    )

    # Reference lines at μ = 1/3, 3 (admissibility transitions) and μ = 1
    ax.axvline(x=1 / 3, color="white", ls="--", lw=1.3, alpha=0.9)
    ax.axvline(x=3, color="white", ls="--", lw=1.3, alpha=0.9)
    ax.axvline(x=1, color="#2ecc71", ls=":", lw=1.2, alpha=0.8)

    # Labels on reference lines (inside plot, near top)
    top = np.pi / 2 * 0.94
    ax.text(
        1 / 3 + 0.12,
        top,
        r"$\frac{1}{3}$",
        ha="left",
        fontsize=9,
        color="white",
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.12", facecolor="#444", edgecolor="none", alpha=0.8
        ),
    )
    ax.text(
        3 + 0.12,
        top,
        r"$3$",
        ha="left",
        fontsize=9,
        color="white",
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.12", facecolor="#444", edgecolor="none", alpha=0.8
        ),
    )
    ax.text(
        1 + 0.12,
        top,
        r"$1$",
        ha="left",
        fontsize=9,
        color="#2ecc71",
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.12", facecolor="#1a1a2e", edgecolor="none", alpha=0.6
        ),
    )

    # Annotate the fully-admissible corridor
    ax.text(
        1.65,
        np.pi / 4,
        "All $\\theta$\nadmissible",
        fontsize=9,
        ha="center",
        color="white",
        fontweight="bold",
        style="italic",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="#333", edgecolor="none", alpha=0.45
        ),
    )

    # Mark the equilateral (theta=0) special point
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

    # Annotate restricted windows
    ax.text(
        5.0,
        np.pi / 3,
        "Restricted\nangular\nwindows",
        fontsize=8,
        ha="center",
        color="#555",
        style="italic",
    )

    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\theta$", rotation=0, labelpad=10)
    ax.set_xlim(0, 6)
    ax.set_ylim(0, np.pi / 2)
    ax.set_box_aspect(1)  # square axes box, matching Mathematica's default
    # x-ticks matching Mathematica: 0, 1/3, 3, 6
    ax.set_xticks([0, 1 / 3, 3, 6])
    ax.set_xticklabels([r"$0$", r"$\frac{1}{3}$", r"$3$", r"$6$"])
    # y-ticks matching Mathematica: 0, π/6, π/3, π/2
    ax.set_yticks([0, np.pi / 6, np.pi / 3, np.pi / 2])
    ax.set_yticklabels(
        [r"$0$", r"$\frac{\pi}{6}$", r"$\frac{\pi}{3}$", r"$\frac{\pi}{2}$"]
    )

    # Legend
    legend_elements = [
        mpatches.Patch(
            color="#e74c3c",
            label=r"$\alpha_{12}$ binding (apex $\leftrightarrow$ base-left)",
        ),
        mpatches.Patch(
            color="#e67e22",
            label=r"$\alpha_{13}$ binding (apex $\leftrightarrow$ base-right)",
        ),
        mpatches.Patch(color="#7f8c8d", label=r"$\alpha_{23}$ binding (base edge)"),
        mpatches.Patch(
            color="#ececec", label="Inadmissible ($\\exists\\,\\alpha_{ij} < 0$)"
        ),
        Line2D([0], [0], color="#1a1a2e", lw=2, label="Admissibility boundary"),
    ]

    ax.legend(
        handles=legend_elements,
        loc="upper right",
        fontsize=8,
        framealpha=0.95,
        edgecolor="#ccc",
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "figure_1_admissibility_binding.pdf")
    fig.savefig(OUTPUT_DIR / "figure_1_admissibility_binding.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
