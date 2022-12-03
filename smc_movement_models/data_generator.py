"""Module pour générer des données similaire à celles de l'article."""

import smc_movement_models.graph_values as gv
import numpy as np


def generate_ts(min=0, max=24.01, step=0.05, num=3000):
    if num is None:
        ts = np.arange(min, max, step)
    else:
        ts = np.linspace(min, max, num)
    return ts


def interpolate_graphs(ts, a1=gv.a1, a2=gv.a2, se=gv.se, s0=gv.s0):
    a1s = np.interp(x=ts, xp=list(a1.keys()), fp=list(a1.values()))
    a2s = np.interp(x=ts, xp=list(a2.keys()), fp=list(a2.values()))
    ses = np.interp(x=ts, xp=list(se.keys()), fp=list(se.values()))
    s0s = np.interp(x=ts, xp=list(s0.keys()), fp=list(s0.values()))

    return a1s, a2s, ses, s0s


def generate_x_tilde(
    ts,
    ses,
    s0s,
    a10=gv.a10,
    a20=gv.a20,
    z0=gv.z0,
    z1=gv.z1,
    zeta0=gv.zeta0,
    c1=gv.c1,
    c2=gv.c2,
    delta=gv.delta,
    sv=gv.sv,
    seed=123,
):
    np.random.seed(seed)

    n = len(ts)
    epst = c1 * np.random.normal(loc=0, scale=np.sqrt(ses)) + c2 * np.random.normal(
        loc=0, scale=delta * np.sqrt(ses)
    )
    nu1t = np.random.normal(loc=0, scale=np.sqrt(sv), size=n)
    nu2t = np.random.normal(loc=0, scale=np.sqrt(sv), size=n)
    nt = np.vstack((epst, np.zeros_like(ts), nu1t, nu2t)).T

    zeta1 = z0 + epst[0]
    a11 = a10 + nu1t[1]
    a21 = a20 + nu2t[1]

    xt = [[z0 + epst[0], zeta0, a10, a20], [z1 + epst[1], zeta1, a11, a21]]

    for t in range(2, n):
        Gt = np.array([[xt[-1][2], xt[-1][3], 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        x = Gt @ np.array(xt[-1]) + nt[t]

        xt.append(list(x))

    xt = np.array(xt)

    et = np.random.normal(loc=0, scale=np.sqrt(s0s))
    yt = xt[:, 0] + et

    return yt, et, xt, nt
