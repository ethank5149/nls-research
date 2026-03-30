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
        2, 1, figsize=(9, 6), height_ratios=[3, 1], gridspec_kw={"hspace": 0.25}
    )

    # Main plot
    ax1.fill_between(x, 0, u_sum, color="#4a90d9", alpha=0.12)  # type: ignore[arg-type]
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

    # Peak positions
    for zi, lab, col in [
        (z1, r"$z_1$", "#e74c3c"),
        (z2, r"$z_2$", "#27ae60"),
        (z3, r"$z_3$", "#f39c12"),
    ]:
        ax1.axvline(zi, color=col, ls=":", lw=0.6, alpha=0.4)
        ax1.text(zi, -0.1, lab, ha="center", fontsize=10, color=col)

    # Minimum inter-peak distance
    ax1.annotate(
        "",
        xy=(z2, 0.5),
        xytext=(z1, 0.5),
        arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9),
    )
    ax1.text(
        (z1 + z2) / 2,
        0.55,
        r"$m_R = \min_{i \neq j} |z_i - z_j|$",
        ha="center",
        fontsize=9,
        color="#555",
    )

    # Profile equation annotation
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

    ax1.set_ylabel(r"$u_E(x)$")
    ax1.set_xlim(-14, 14)
    ax1.set_ylim(-0.15, 1.25)
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax1.set_title(
        r"Multi-Peak Decomposition: $u_E(x) \approx \sum_{i=1}^{M} u_R(x - z_i(R))$",
        fontsize=12,
    )

    # Remainder
    ax2.plot(x, h_R, color="#8e44ad", lw=1.2)
    ax2.axhline(0, color="k", lw=0.3)
    ax2.fill_between(x, 0, h_R, color="#8e44ad", alpha=0.12)  # type: ignore[arg-type]
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"$h_R(x)$")
    ax2.set_xlim(-14, 14)
    ax2.set_title(
        r"Remainder $h_R$: $\|h_R\|_{H^1} \to 0$ as $E \to \infty$ "
        r"(solved via IFT, Lemma 2)",
        fontsize=10,
    )

    fig.savefig(OUTPUT_DIR / "figure_4_profile_decomposition.pdf")
    fig.savefig(OUTPUT_DIR / "figure_4_profile_decomposition.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
