"""Module pour créer des figures."""

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import seaborn as sns
import numpy as np

plt.style.use("seaborn-v0_8-paper")


def plot_acvf(acvfs, window_times, figsize=(10, 5)):
    x_lims = window_times[0], window_times[-1]
    x_lims = mdates.date2num(x_lims)
    y_lims = [60, 0]

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    im = ax.imshow(
        acvfs, extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]], aspect="auto", cmap="Spectral_r"
    )
    ax.set(ylabel="Lag (s)", xlabel="Time of day")
    ax.xaxis_date()
    date_format = mdates.DateFormatter("%H:%M:%S")
    ax.xaxis.set_major_formatter(date_format)
    fig.autofmt_xdate()

    fig.colorbar(im)

    return fig


def plot_variance_estimates(sigma_eds, sigma_ods, window_times, figsize=(7, 3.5)):
    res_estimates = pd.DataFrame(
        {
            "Dtime": window_times,
            "sigma_ed": sigma_eds,
            "sigma_od": sigma_ods,
        }
    )

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    res_estimates.plot(
        x="Dtime",
        y="sigma_od",
        color="tab:blue",
        marker="^",
        markersize=5,
        label=r"$\hat\sigma_o^2$ = observation error variance",
        ax=ax,
    )
    res_estimates.plot(
        x="Dtime",
        y="sigma_ed",
        color="tab:red",
        marker="o",
        markersize=5,
        label=r"$\hat\sigma_{\varepsilon}^2$ = system error variance",
        ax=ax,
    )
    ax.set(xlabel="Time of day", ylabel="Variance")

    return fig


def plot_depth(data, figsize=(7, 3.5)):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    data.plot(x="Dtime", y="Depth", ax=ax)
    ax.set(xlabel="Time of day", ylabel="Depth (m)")
    return fig


def plot_velocity(data, figsize=(7, 3.5)):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    data.plot(x="Dtime", y="Velocity", ax=ax)
    ax.set(xlabel="Time of day", ylabel="Vertical velocity (m/s)")
    return fig


def plot_a1_a2(time_results, figsize=(10, 10)):
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=figsize)
    time_results.plot(
        x="Dtime", y="a1s", ax=axs[0], marker="o", color="moccasin", label=r"$\hat{a}_1$"
    )
    time_results.plot(
        x="Dtime", y="rolling_a1s", ax=axs[0], color="orange", label="Rolling average"
    )
    axs[0].set(ylabel=r"$\hat{a}_1$", xlabel="Time of day")
    axs[0].legend(loc="upper left")
    time_results.plot(
        x="Dtime", y="a2s", ax=axs[1], marker="^", color="lightblue", label=r"$\hat{a}_2$"
    )
    time_results.plot(x="Dtime", y="rolling_a2s", ax=axs[1], color="blue", label="Rolling average")
    axs[1].set(ylabel=r"$\hat{a}_2$", xlabel="Time of day")
    fig.tight_layout()
    return fig


def plot_a1_a2_single_plot(time_results, figsize=(7, 3.5)):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    time_results.plot(x="Dtime", y="a1s", ax=ax, marker="o", color="moccasin", label=r"$\hat{a}_1$")
    time_results.plot(x="Dtime", y="rolling_a1s", ax=ax, color="orange", label="Rolling average")
    time_results.plot(
        x="Dtime", y="a2s", ax=ax, marker="^", color="lightblue", label=r"$\hat{a}_2$"
    )
    time_results.plot(x="Dtime", y="rolling_a2s", ax=ax, color="blue", label="Rolling average")
    ax.set(ylabel="Variances", xlabel="Time of day")
    # ax.legend(loc="upper left")
    return fig


def plot_a1_a2_var(final_res, figsize=(7, 3.5)):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    sns.lineplot(data=final_res, x="Dtime", y="a1s", ax=ax, label=r"$\hat{a}_1$")
    sns.lineplot(data=final_res, x="Dtime", y="a2s", ax=ax, label=r"$\hat{a}_2$")
    # fig.legend()
    ax.set(xlabel="Time of day", ylabel="Variances")
    x_ids = np.linspace(0, 108, 9)
    ax.set_xticklabels(
        labels=final_res["Dtime"].iloc[x_ids].dt.strftime("%H:%M"), rotation=45, ha="right"
    )
    return fig
