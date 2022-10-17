import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from lib.series import convert_to_series

plt.rc("text", usetex=True)


series = convert_to_series(0, 4.2, 0.01, [1, 0.5, 1], [1, 1.5, 1])

with PdfPages("figs/explanation.pdf") as pdfFile:
    fig = plt.figure(figsize=(3, 2))

    ax = plt.gca()
    ax.grid(visible=True, linestyle="--", dashes=(15, 15), c="#cccccc", lw=0.3)
    ax.get_yaxis().set_tick_params(direction="out")
    ax.get_xaxis().set_tick_params(direction="out")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$I(t)$")
    ax.set_yticks([0, 1])
    ax.set_yticklabels([r"$0$", r"$a$"])
    ax.set_xticks([0, 1, 3.5])
    ax.set_xticklabels([r"$t=0$", r"$t_1$", r"$t_2$"])

    plt.plot(series[:, 0], series[:, 1], "r-")
    plt.plot([4.25, 5], [0, 0], "r:")

    ax.annotate(r"$\tau_0$", [0.5, 0.5], ha="center", va="center", xycoords="data")
    ax.annotate(
        "",
        xy=(0, 0.4),
        xycoords="data",
        xytext=(1, 0.4),
        textcoords="data",
        arrowprops=dict(arrowstyle="<->"),
    )
    ax.annotate(r"$\theta_1$", [1.5, 0.5], ha="center", va="center", xycoords="data")
    ax.annotate(
        "",
        xy=(1, 0.6),
        xycoords="data",
        xytext=(2, 0.6),
        textcoords="data",
        arrowprops=dict(arrowstyle="<->"),
    )
    ax.annotate(r"$\tau_1$", [2.75, 0.5], ha="center", va="center", xycoords="data")
    ax.annotate(
        "",
        xy=(2, 0.4),
        xycoords="data",
        xytext=(3.5, 0.4),
        textcoords="data",
        arrowprops=dict(arrowstyle="<->"),
    )
    ax.annotate(r"$\theta_2$", [3.75, 0.5], ha="center", va="center", xycoords="data")
    ax.annotate(
        "",
        xy=(3.5, 0.6),
        xycoords="data",
        xytext=(4, 0.6),
        textcoords="data",
        arrowprops=dict(arrowstyle="<->"),
    )

    plt.tight_layout(w_pad=1, h_pad=1, pad=0.5)

    pdfFile.savefig(fig)
    plt.close()
