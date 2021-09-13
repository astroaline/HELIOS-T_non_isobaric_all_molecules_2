import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm

fig, ax = pl.subplots()

mu_rstar = 1.57
sigma_rstar =  0.07
x_gauss_rstar = np.linspace(mu_rstar - 3*sigma_rstar, mu_rstar + 3*sigma_rstar, 100)
pl.plot(x_gauss_rstar, 90*norm.pdf(x_gauss_rstar, mu_rstar, sigma_rstar), color='r', linestyle='-', linewidth=4)
pl.xlim((1.5,1.65))
pl.xlabel('R$_{star}$')


#mu_g = 977
#sigma_g = 67
#x_gauss_g = np.linspace(mu_g - 3*sigma_g, mu_g + 3*sigma_g, 100)
#pl.plot(x_gauss_g, 90*norm.pdf(x_gauss_g, mu_g, sigma_g), color='g', linestyle='-', linewidth=4)
#pl.xlim((900,1050))
#pl.xlabel('g')

pl.savefig("rstar_gauss.png", transparent=False, bbox_inches='tight')
