import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

plt.rc("text", usetex=True)


with PdfPages("figs/poiss-poiss.pdf") as pdfFile:
    data = 10 ** np.loadtxt("data/poiss10.poiss10.seed7267.psd.csv", delimiter=",")

    fig = plt.figure(figsize=(3, 2))
    plt.loglog()

    ax = plt.gca()
    ax.grid(visible=True, linestyle="--", dashes=(15, 15), c="#cccccc", lw=0.3)
    ax.get_yaxis().set_tick_params(direction="out")
    ax.get_xaxis().set_tick_params(direction="out")
    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$S(f)$")

    plt.plot(data[:, 0], data[:, 1], "r-", lw=4)
    plt.plot(data[:, 0], data[:, 2], "k-")

    plt.tight_layout(w_pad=1, h_pad=1, pad=0.5)

    pdfFile.savefig(fig)
    plt.close()
