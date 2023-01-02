"""Module qui implemente le modèle présenté dans le papier principal.

Examples:
    Exemple d'implementation avec `particles`_:

    ```python
        import particles
        import particles.state_space_models as ssm
        import particles.distributions as dists

        class ToySSM(ssm.StateSpaceModel):
            def PX0(self):  # Distribution of X_0
                return dists.Normal()  # X_0 ~ N(0, 1)
            def PX(self, t, xp):  # Distribution of X_t given X_{t-1}
                return dists.Normal(loc=xp)  # X_t ~ N( X_{t-1}, 1)
            def PY(self, t, xp, x):  # Distribution of Y_t given X_t (and X_{t-1})
                return dists.Normal(loc=x, scale=self.sigma)  # Y_t ~ N(X_t, sigma^2)
    ```

.. _particles:
    https://github.com/nchopin/particles
"""

from particles.distributions import ProbDist, Categorical

import scipy as sp
import numpy as np

from particles import SMC
from particles.distributions import Normal, IndepProd
from particles.state_space_models import Bootstrap, StateSpaceModel

from tqdm import tqdm
import pandas as pd
from statsmodels.tsa.stattools import acovf


class Mixture(ProbDist):  # Copied from the last commit of particles (not yet released).
    """Mixture distributions.
    Parameters
    ----------
    pk:  array-like of shape (k,) or (N, k)
        component probabilities (must sum to one)
    *components: ProbDist objects
        the k component distributions
    Example:
        mix = Mixture([0.6, 0.4], Normal(loc=3.), Normal(loc=-3.))
    """

    def __init__(self, pk, *components):
        self.pk = np.atleast_1d(pk)
        self.k = self.pk.shape[-1]  # Fixed typo
        if len(components) != self.k:
            raise ValueError("Size of pk and nr of components should match")
        self.components = components

    def logpdf(self, x):
        lpks = [np.log(self.pk[..., i]) + cd.logpdf(x) for i, cd in enumerate(self.components)]
        return sp.special.logsumexp(np.column_stack(tuple(lpks)), axis=-1)

    def rvs(self, size=None):
        k = Categorical(p=self.pk).rvs(size=size)
        xk = [cd.rvs(size=size) for cd in self.components]
        # sub-optimal, we sample N x k values
        return np.choose(k, xk)


class MarineSSM(StateSpaceModel):
    """State space model for the seal model."""

    default_params = {
        "z0": 0.0,  # Initial value of z_0
        "z1": 0.0,  # Initial value of z_1
        "a1": 0.0,  # Initial value of a_1
        "a2": 0.0,  # Initial value of a_2
        "sigma_o": 0.1,  # observation error (y) <-> e
        "sigma_v": 0.1,  # Disturbance term (a1, a2) <-> nu
        "sigma_e": 0.1,  # System noise (z) <-> epsilon
        "c1": 0.9,  # Mixture options for the system noise
        "c2": 0.1,
        "delta": 10,
    }

    def PX0(self):  # dist of X_0
        zt = Mixture(  # z_1 (mixture)
            [self.c1, self.c2],
            Normal(loc=self.z1, scale=self.sigma_e),
            Normal(loc=self.z1, scale=self.delta * self.sigma_e),
        )
        ztm = Mixture(  # z_0 (mixture)
            [self.c1, self.c2],
            Normal(loc=self.z0, scale=self.sigma_e),
            Normal(loc=self.z0, scale=self.delta * self.sigma_e),
        )
        a1t = Normal(loc=self.a1, scale=self.sigma_v)  # a1_t
        a2t = Normal(loc=self.a2, scale=self.sigma_v)  # a2_t
        return IndepProd(zt, ztm, a1t, a2t)

    def PX(self, t, xp):  # dist of X_t at time t, given X_{t-1}
        loc = xp[:, 2] * xp[:, 0] + xp[:, 3] * xp[:, 1]
        zt = Mixture(  # z_t (mixture)
            [self.c1, self.c2],
            Normal(loc=loc, scale=self.sigma_e),
            Normal(loc=loc, scale=self.delta * self.sigma_e),
        )
        ztm = Normal(loc=xp[:, 0], scale=0.0)  # z_{t-1} (0 variance)
        a1t = Normal(loc=xp[:, 2], scale=self.sigma_v)  # a1_t
        a2t = Normal(loc=xp[:, 3], scale=self.sigma_v)  # a2_t
        return IndepProd(zt, ztm, a1t, a2t)

    def PY(self, t, xp, x):  # dist of Y_t at time t, given X_t and X_{t-1}
        return Normal(loc=x[:, 0], scale=self.sigma_o)


def get_windows(day, overlap_size=156):
    """Generate overlapping windows of size 2 * overlap_size - 1 from a day of data."""
    num_windows = int((len(day) - overlap_size) / overlap_size)
    windows = []
    window_times = []
    for i in range(num_windows):
        window = day[i * overlap_size : (i + 2) * overlap_size - 1]["Velocity"].values
        window_time = day["Dtime"].iloc[(i + 1) * overlap_size]
        windows.append(window)
        window_times.append(window_time)
    return windows, window_times


def estimate_variances(window, nlag=60, n_estimates=8):
    """Estimate the variances of the observation, system noise for a window."""
    acvf = acovf(window, nlag=nlag)
    x = np.arange(n_estimates + 1)
    y = acvf[: n_estimates + 1]
    coeffs = np.polyfit(x[1:], y[1:], 2)
    regr = np.poly1d(coeffs)
    mse = np.mean((regr(x) - y) ** 2)

    sigma_zd = np.max([coeffs[-1], 1e-4])
    sigma_od = np.max([acvf[0] - sigma_zd, 1e-4])
    sigma_ed = (1 - (acvf[1] / acvf[0]) ** 2) * sigma_zd

    return sigma_od, sigma_ed, sigma_zd, mse, acvf


def estimate_variance_on_all_windows(
    windows,
    nlag=12,
    n_estimates=8,
):
    """Estimate the variances of the observation, system noise for all
    windows."""
    sigma_ods = []
    sigma_eds = []
    sigma_zds = []
    mses = []
    acvfs = []
    for window in windows:
        sigma_od, sigma_ed, sigma_zd, mse, acvf = estimate_variances(
            window, nlag=nlag, n_estimates=n_estimates
        )
        sigma_ods.append(sigma_od)
        sigma_eds.append(sigma_ed)
        sigma_zds.append(sigma_zd)
        mses.append(mse)
        acvfs.append(acvf)

    sigma_ods = np.array(sigma_ods)
    sigma_eds = np.array(sigma_eds)
    sigma_zds = np.array(sigma_zds)
    mses = np.array(mses)
    acvfs = np.array(acvfs).T

    return sigma_ods, sigma_eds, sigma_zds, mses, acvfs


def run_smc(window, N=500, **kwargs):
    """Run the SMC on a window (one iteration of the algorithm)."""
    my_ssm_model = MarineSSM(z0=window[0], z1=window[1], **kwargs)
    my_fk_model = Bootstrap(ssm=my_ssm_model, data=window)
    my_alg = SMC(fk=my_fk_model, N=N, store_history=True)
    my_alg.run()
    return my_alg


def estimate_a1_a2_on_window(
    window,
    N=500,
    M=10,
    initial_a1=0,
    initial_a2=0,
    initial_sigma_v=0.1,
    sigma_o=0.1,
    sigma_e=0.1,
    alpha=0.5,
    c1=0.9,
    c2=0.1,
    delta=10,
):
    """Run the algorithm on a window."""
    final_a1s = []
    final_a2s = []
    for m in range(M):
        a1 = final_a1s[-1].mean() if m > 0 else initial_a1
        a2 = final_a2s[-1].mean() if m > 0 else initial_a2
        alg = run_smc(
            window=window,
            N=N,
            sigma_v=initial_sigma_v,
            a1=a1,
            a2=a2,
            sigma_o=sigma_o,
            sigma_e=sigma_e,
            c1=c1,
            c2=c2,
            delta=delta,
        )
        a1s_tmp = np.array(alg.hist.X)[:, :, 2]  # shape (T, N)
        a2s_tmp = np.array(alg.hist.X)[:, :, 3]
        wgts = np.array([w.W for w in alg.hist.wgts])  # shape (T, N)
        final_a1s.append(np.average(a1s_tmp, weights=wgts, axis=1))  # shape (T)
        final_a2s.append(np.average(a2s_tmp, weights=wgts, axis=1))
        initial_sigma_v *= alpha

    final_a1s = np.array(final_a1s)
    final_a2s = np.array(final_a2s)

    return final_a1s, final_a2s


def estimate_a1_a2_on_all_windows(
    windows,
    window_times,
    sigma_os,
    sigma_es,
    N=500,
    M=10,
    initial_a1=0,
    initial_a2=0,
    initial_sigma_v=0.1,
    alpha=0.5,
    c1=0.9,
    c2=0.1,
    delta=10,
):
    """Run the algorithm on all windows."""
    time_a1s = []
    time_a2s = []

    for window, sigma_o, sigma_e in tqdm(zip(windows, sigma_os, sigma_es), total=len(windows)):
        final_a1s, final_a2s = estimate_a1_a2_on_window(
            window=window,
            N=N,
            M=M,
            initial_a1=initial_a1,
            initial_a2=initial_a2,
            initial_sigma_v=initial_sigma_v,
            sigma_o=sigma_o,
            sigma_e=sigma_e,
            alpha=alpha,
            c1=c1,
            c2=c2,
            delta=delta,
        )
        time_a1s.append(final_a1s)
        time_a2s.append(final_a2s)

    time_a1s = np.array(time_a1s)  # shape (n_windows, M, T)
    time_a2s = np.array(time_a2s)
    time_results = pd.DataFrame(  # shape (n_windows, 5)
        {
            "Dtime": window_times,
            "a1s": time_a1s[:, -1].mean(axis=1),
            "a2s": time_a2s[:, -1].mean(axis=1),
        }
    )
    time_results["rolling_a1s"] = time_results["a1s"].rolling(window=5, center=True).mean()
    time_results["rolling_a2s"] = time_results["a2s"].rolling(window=5, center=True).mean()

    return time_results


def rerun_algo(
    windows,
    window_times,
    sigma_os,
    sigma_es,
    nb_runs=20,
    **kwargs,
):
    results_var = []

    for i in tqdm(range(nb_runs)):
        time_results = estimate_a1_a2_on_all_windows(
            windows=windows,
            window_times=window_times,
            sigma_os=sigma_os,
            sigma_es=sigma_es,
            **kwargs,
        )
        results_var.append(time_results)

    final_res = pd.DataFrame()
    for i, res in enumerate(results_var):
        res["run"] = i
        final_res = pd.concat([final_res, res])

    return final_res
