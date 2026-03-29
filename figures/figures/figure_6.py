import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
    print("Creating Figure 6 - Interaction Functional Decay...", end=" ", flush=True)
    d = np.linspace(0.5, 12, 500)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    R_vals = [0.01, 0.1, 0.3, 0.5, 0.8]
    V0 = 1.0
    colors = cm.get_cmap("plasma")(np.linspace(0.15, 0.85, len(R_vals)))

    for R, col in zip(R_vals, colors):
        kappa = np.sqrt(1 + R**2 * V0)
        Q = np.exp(-d * kappa) * d ** (-0.5)
        ax1.semilogy(d, Q, color=col, lw=1.8, label=rf"$R = {R}$")

    ax1.set_xlabel(r"Inter-peak distance $d = |z_k - z_i|$")
    ax1.set_ylabel(r"$Q_R(d)$ (log scale)")
    ax1.legend(fontsize=9, title=r"$R = E^{-1/2}$", title_fontsize=9)
    ax1.set_xlim(0.5, 12)
    ax1.set_ylim(1e-8, 1)

    # Annotate decay rate
    ax1.annotate(
        r"$Q_R(d) \sim C \cdot d^{-(n-1)/2} \cdot e^{-d\sqrt{1 + R^2 V(x_0)}}$",
        xy=(3, 1e-2),
        fontsize=9,
        color="#333",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )

    # Shade negligible region
    ax1.axhspan(0, 1e-6, color="#f0f0f0", alpha=0.6, zorder=0)
    ax1.text(
        8.5,
        3e-7,
        "Negligible interaction\n(justifies $N_i$ truncation)",
        fontsize=8,
        color="#7f8c8d",
        ha="center",
        style="italic",
    )

    # Right panel: α_ki vs distance ratio
    d_ratio = np.linspace(1, 3.5, 500)
    for R, col in zip(R_vals, colors):
        kappa = np.sqrt(1 + R**2 * V0)
        m_R = 3.0
        alpha = np.exp(-m_R * kappa * (d_ratio - 1)) * d_ratio ** (-0.5)
        ax2.plot(d_ratio, alpha, color=col, lw=1.8, label=rf"$R = {R}$")

    ax2.axhline(1, color="k", ls=":", lw=0.5)
    ax2.axhline(0, color="k", lw=0.3)
    ax2.set_xlabel(r"Distance ratio $|z_k - z_i| \, / \, m_R$")
    ax2.set_ylabel(r"$\alpha_{ki} = Q_R(|z_k - z_i|) \, / \, Q_R(m_R)$")
    ax2.set_xlim(1, 3.5)
    ax2.set_ylim(-0.05, 1.15)

    # Annotate regimes
    ax2.annotate(
        "Nearest neighbors\n" r"$|z_k - z_i| = m_R$" "\n" r"$\alpha_{ki} = 1$",
        xy=(1.02, 0.98),
        xytext=(1.4, 0.65),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="#555", lw=0.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#d5f5e3", edgecolor="#82e0aa"),
    )
    ax2.annotate(
        "Next-nearest neighbors\n" r"$\alpha_{ki} \ll 1$" "\n(dropped from sum)",
        xy=(2.0, 0.03),
        xytext=(2.3, 0.4),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="#555", lw=0.8),
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#fadbd8", edgecolor="#f1948a"),
    )

    # Annotation: definition of α
    ax2.text(
        0.02,
        0.02,
        r"$0 \leq \alpha_{ki} \leq 1$ by construction; "
        "\n"
        r"$\alpha_{ki} = 1$ iff $k \in N_i$",
        transform=ax2.transAxes,
        fontsize=8,
        va="bottom",
        color="#555",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="#fef9e7", edgecolor="#ddd"),
    )
    ax2.legend(fontsize=8, loc="center right")

    fig.tight_layout()
    fig.savefig(f"figures/output/figure_6_interaction_decay.pdf")
    fig.savefig(f"figures/output/figure_6_interaction_decay.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
