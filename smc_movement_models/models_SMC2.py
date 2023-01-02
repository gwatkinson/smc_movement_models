"""Module qui implemente le modèle SMC2."""

from particles.distributions import StructDist

import numpy as np

from particles import SMC
from particles.smc_samplers import SMC2
from particles.distributions import Normal, IndepProd
from particles.state_space_models import StateSpaceModel
from smc_movement_models.models import Mixture

fenetre = [0, 0]


class MarineSSM_SMC2(StateSpaceModel):
    """Classe qui implemente le state space model pour le SMC2."""

    global fenetre
    default_params = {
        "z0": fenetre[0],  # Initial value of z_0
        "z1": fenetre[1],  # Initial value of z_1
        "a1": 0.0,  # Value of a_1
        "a2": 0.0,  # Value of a_2
        "sigma_o": 0.5,  # observation error (y) <-> e
        "sigma_e": 0.5,  # System noise (z) <-> epsilon
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


def run_smc2(window, N=10, sigma_1=0.1, sigma_2=0.1, init_Nx=5, ar_to_increase_Nx=0.1):
    """Run l'algorithme SMC2 sur une fenetre de données."""
    global fenetre
    fenetre = window
    my_ssm_model = MarineSSM_SMC2
    my_prior = StructDist({"a1": Normal(scale=sigma_1), "a2": Normal(scale=sigma_2)})
    my_fk_model = SMC2(
        ssm_cls=my_ssm_model,
        prior=my_prior,
        data=window,
        init_Nx=init_Nx,
        ar_to_increase_Nx=ar_to_increase_Nx,
    )
    my_alg = SMC(fk=my_fk_model, N=N, store_history=True)
    my_alg.run()
    return my_alg


def estimate_a1_a2_on_window_SMC2(
    window, N=10, sigma_1=0.1, sigma_2=0.1, init_Nx=5, ar_to_increase_Nx=0.1
):
    """Run l'algorithme SMC2 sur toute les fenetres de données."""
    alg = run_smc2(
        window=window,
        N=N,
        sigma_1=sigma_1,
        sigma_2=sigma_2,
        init_Nx=init_Nx,
        ar_to_increase_Nx=ar_to_increase_Nx,
    )
    wgts = np.array([w.W for w in alg.hist.wgts])
    a_temp = np.array([x.theta for x in alg.hist.X])
    a1s_temp = a_temp["a1"]
    a2s_temp = a_temp["a2"]
    return np.average(a1s_temp, weights=wgts, axis=1), np.average(a2s_temp, weights=wgts, axis=1)
