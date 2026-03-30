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
    """
    Illustrates Claim 2: Rz_i(R) → x₀.
    """
    print(
        "Creating Figure 5 - Concentration Points vs Peak Positions...",
        end=" ",
        flush=True,
    )
    fig = plt.figure(figsize=(14, 8))

    # Left: 2x2 grid of profiles at different E
    # Right: convergence plot of peak positions
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.2], hspace=0.35, wspace=0.35)

    ax_profiles = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    ax_conv = fig.add_subplot(gs[:, 2])

    x = np.linspace(-10, 10, 1000)
    x0 = 0.0  # critical point

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

    for ax, cfg in zip(ax_profiles, configs):
        psi = np.zeros_like(x)
        for z, _, w in cfg["peaks"]:
            psi += sech(x, z, w)

        ax.fill_between(x, 0, psi, color="#4a90d9", alpha=0.15)  # type: ignore[arg-type]
        ax.plot(x, psi, color="#1a3a6e", lw=2.0)

        # Mark x₀
        ax.axvline(x0, color="#27ae60", ls="--", lw=1.2, alpha=0.7)
        ax.text(x0 + 0.3, max(psi) * 0.92, r"$x_0$", fontsize=10, color="#27ae60")

        # Mark peak positions z_i
        for i, (z, _, w) in enumerate(cfg["peaks"]):
            ax.plot(z, sech(z, z, w), "v", color="#e74c3c", ms=8, zorder=5)
            ax.text(z, -0.12, rf"$z_{i + 1}$", ha="center", fontsize=9, color="#e74c3c")

        # Mark Rz_i values (the rescaled positions converging to x₀)
        R = cfg["R"]
        for i, (z, _, w) in enumerate(cfg["peaks"]):
            Rz = R * z
            # Small marker on x-axis for Rz_i
            ax.plot(Rz, -0.05, "s", color="#8e44ad", ms=5, zorder=5, clip_on=False)

        ax.axhline(0, color="k", lw=0.3)
        ax.set_xlim(-8, 8)
        ax.set_ylim(-0.2, 1.3)
        ax.set_title(cfg["title"], fontsize=10)
        ax.set_xlabel(r"$x$" if ax in ax_profiles[2:] else "")
        ax.set_ylabel(r"$|\psi_E(x)|$")

        # Info box
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

    # Legend for profile plots
    ax_profiles[0].plot(
        [], [], "v", color="#e74c3c", ms=7, label=r"Peak centers $z_i(E)$"
    )
    ax_profiles[0].plot([], [], "s", color="#8e44ad", ms=5, label=r"Rescaled $Rz_i$")
    ax_profiles[0].plot(
        [], [], "--", color="#27ae60", lw=1.2, label=r"Critical point $x_0$"
    )
    ax_profiles[0].legend(fontsize=7, loc="upper left", framealpha=0.9)

    # ── Right panel: convergence of Rz_i → x₀ ──
    E_vals = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    R_vals = E_vals ** (-0.5)

    # Simulate M=2 peak positions: z_i ≈ ±C/R as R→0 so Rz_i → ±C·R/R... no
    # Actually: Rz_i → x₀ as R→0, so z_i ~ x₀/R + small correction
    # Let's model: Rz_1(R) = x₀ + δ(R) where δ(R) → 0
    # For two peaks: Rz_± ≈ x₀ ± c·R for some constant c

    np.random.seed(42)
    c1, c2 = 3.5, -3.5  # asymptotic peak spread
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

    # Shade convergence corridor
    ax_conv.fill_between(E_vals, x0 - 0.05, x0 + 0.05, color="#27ae60", alpha=0.1)

    ax_conv.set_xlabel(r"Energy $E$", fontsize=12)
    ax_conv.set_ylabel(r"Rescaled position $Rz_i(E) = E^{-1/2}z_i(E)$", fontsize=11)
    ax_conv.set_title(r"$Rz_i(E) \to x_0$ as $E \to \infty$" "\n(Claim 2)", fontsize=12)
    ax_conv.legend(fontsize=9, loc="upper right")
    ax_conv.set_ylim(-2.5, 2.5)
    ax_conv.grid(True, alpha=0.2, which="both")

    # Annotation
    ax_conv.annotate(
        r"$|Rz_i - x_0| = O(E^{-1/2})$",
        xy=(100, 0.4),
        fontsize=10,
        color="#555",
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fef9e7", edgecolor="#ddd"),
    )

    fig.savefig(OUTPUT_DIR / "figure_5_concentration_vs_peaks.pdf")
    fig.savefig(OUTPUT_DIR / "figure_5_concentration_vs_peaks.png")
    plt.close(fig)
    print("Done!")


if __name__ == "__main__":
    create()
