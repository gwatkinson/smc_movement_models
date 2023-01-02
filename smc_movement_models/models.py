"""Module qui implemente les modèles.

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

from particles.distributions import ProbDist, Categorical, StructDist

import scipy as sp
import numpy as np

from particles import SMC
from particles.smc_samplers import SMC2
from particles.distributions import Normal, IndepProd
from particles.state_space_models import Bootstrap, StateSpaceModel


class Mixture(ProbDist):
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
        self.k = self.pk.shape[-1]
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


class MarineSSM_SMC2(StateSpaceModel):
    default_params = {
        "z0": 0.0,  # Initial value of z_0
        "z1": 0.0,  # Initial value of z_1
        "a1": 0.0,  # Value of a_1
        "a2": 0.0,  # Value of a_2
        "sigma_o": 0.1,  # observation error (y) <-> e
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
        return IndepProd(zt, ztm)

    def PX(self, t, xp):  # dist of X_t at time t, given X_{t-1}
        loc = self.a1 * xp[:, 0] + self.a2 * xp[:, 1]
        zt = Mixture(  # z_t (mixture)
            [self.c1, self.c2],
            Normal(loc=loc, scale=self.sigma_e),
            Normal(loc=loc, scale=self.delta * self.sigma_e),
        )
        ztm = Normal(loc=xp[:, 0], scale=0.0)  # z_{t-1} (0 variance)
        return IndepProd(zt, ztm)

    def PY(self, t, xp, x):  # dist of Y_t at time t, given X_t and X_{t-1}
        return Normal(loc=x[:, 0], scale=self.sigma_o)


def run_smc(window, N=200, **kwargs):
    my_ssm_model = MarineSSM(z0=window[0], z1=window[1], **kwargs)
    my_fk_model = Bootstrap(ssm=my_ssm_model, data=window)
    my_alg = SMC(fk=my_fk_model, N=N, store_history=True)
    my_alg.run()
    return my_alg

def run_smc2(window, N=200, **kwargs):
    my_ssm_model = MarineSSM_SMC2(z0=window[0], z1=window[1], **kwargs)
    my_prior = StructDist({"prior":Normal(scale=0.1)})
    my_fk_model = SMC2(ssm=my_ssm_model, prior=my_prior, data=window)
    my_alg = SMC(fk=my_fk_model, N=N, store_history=True)
    my_alg.run()
    return my_alg 


def estimate_a1_a2_on_window(window, N=500, M=10, epsilon=0.1, alpha=0.5):
    final_a1s = []
    final_a2s = []
    for m in range(M):
        a1 = final_a1s[-1].mean() if m > 0 else 0
        a2 = final_a2s[-1].mean() if m > 0 else 0
        alg = run_smc(
            window=window,
            sigma_v=epsilon,
            a1=a1,
            a2=a2,
            sigma_o=0.1,
            sigma_e=0.1,
            c1=0.9,
            c2=0.1,
            delta=10,
        )
        a1s_tmp = np.array(alg.hist.X)[:, :, 2]
        a2s_tmp = np.array(alg.hist.X)[:, :, 3]
        wgts = np.array([w.W for w in alg.hist.wgts])
        final_a1s.append(np.average(a1s_tmp, weights=wgts, axis=1))
        final_a2s.append(np.average(a2s_tmp, weights=wgts, axis=1))
        epsilon *= 0.5

    final_a1s = np.array(final_a1s)
    final_a2s = np.array(final_a2s)

    return final_a1s, final_a2s


def estimate_variance():
    pass
