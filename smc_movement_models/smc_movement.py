import pandas as pd
import numpy as np
import particles
import particles.state_space_models as ssm
import particles.distributions as dists
df = pd.read_csv('/home/mbenj/smc_movement_models/data/clean_data.csv')
numpy_df = np.asarray(df)

class Mixture(dists.ProbDist):
    """Base class for Mixture Gaussian distributions.
    """
    def __init__(self, loc1 =0., scale1=1., loc2 = 0, scale2 = 1., ratio = 1/2):
        self.distr1 = dists.Normal(loc = loc1, scale = scale1)
        self.distr2 = dists.Normal(loc = loc2, scale = scale2)
        self.ratio = ratio
    
    def pdf(self, x):
        return self.ratio*self.distr1.pdf(x) + (1-self.ratio)*self.distr2.pdf(x)

    def logpdf(self, x):
        return np.log(self.pdf(x))

    def rvs(self, size=None):
        return self.ratio*self.distr1.rvs(size) + (1 - self.ratio)*self.distr2.rvs(size)
    
    def ppf(self, u):
        return self.ratio*self.distr1.ppf(u) + (1 - self.ratio)*self.distr2.ppf(u)

class ToySSM(ssm.StateSpaceModel):

    default_params = {'c1':.1, 'c2':.9, 'sigma_v':0.1, 'sigma_eps':0.1, 'delta': 10, 'sigma_erreur': 0.1}
    def v_t(self):
        return dists.MvNormal(loc = np.zeros(2), scale = self.sigma_v, cov = np.eye(2))
    def f(self,t, x):
        F = np.zeros((4,4))
        v = self.v_t().rvs(1)
        F[0,:2] = x[0,2] + v[0,0],x[0,3]  + v[0,1]
        F[2:,-1] = 1,1
        return F@x.T, v
    def PX0(self):  # Distribution of X_0 
        return dists.MvNormal(loc = np.zeros(4), cov = np.eye(4))  # X_0 ~ N(0, 1)
    def PX(self, t, xp):  # Distribution of X_t given X_{t-1}
        
        mean, var_v = self.f(t,xp)
        z_t = Mixture(loc1 = mean[0], loc2 = mean[0], scale1 = self.sigma_eps, scale2 = self.delta*self.sigma_eps, ratio = .1)
        ksi_t = dists.Normal(loc = mean[1], scale = 0.0001)
        a1_t = dists.Normal(loc = mean[2] + var_v[0,0], scale = 0.0001)
        a2_t = dists.Normal(loc = mean[3] + var_v[0,0], scale = 0.0001)
        
        return dists.IndepProd(z_t, ksi_t, a1_t, a2_t)  # X_t ~ N( X_{t-1}, 1)
    def PY(self, t, xp, x):  # Distribution of Y_t given X_t (and X_{t-1}) 
        return dists.Normal(loc=x[0,0], scale=self.sigma_erreur)  # Y_t ~ N(X_t, sigma^2)

my_model = ToySSM()
x, y = my_model.simulate(200)
