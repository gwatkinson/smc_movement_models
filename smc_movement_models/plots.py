"""Module pour créer des figures."""

import matplotlib.pyplot as plt
import pandas as pd

import smc_movement_models.graph_values as gv


# Plotting functions
def plot_graph_values(figsize=(10, 10)):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=figsize)

    ax1.plot(
        gv.se.keys(),
        gv.se.values(),
        marker="d",
        label=r"System noise ($\hat{\sigma}_{\varepsilon}^2$)",
    )
    ax1.plot(
        gv.s0.keys(), gv.s0.values(), marker="v", label=r"Observation error ($\hat{\sigma}_{o}^2$)"
    )
    ax1.set(ylabel="Variance")
    ax1.legend()

    ax2.plot(gv.a1.keys(), gv.a1.values())
    ax2.axhline(y=0, ls="--")
    ax2.set(ylabel="$a_1$")

    ax3.plot(gv.a2.keys(), gv.a2.values())
    ax3.axhline(y=0, ls="--")
    ax3.set(ylabel="$a_2$", xlabel="Time of day")

    fig.suptitle("Appoximate results from the paper")
    fig.tight_layout()

    return fig


def plot_real_data(path, figsize=(10, 10)):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize)

    df = pd.read_csv(path, sep="\t")

    df.plot(x="Dtime", y="Velocity", ax=ax1)
    df.plot(x="Dtime", y="Depth", ax=ax2)

    return fig
