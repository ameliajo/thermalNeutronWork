
import numpy as np
from math import factorial as f
import matplotlib.pyplot as plt



alpha = 1.0
betaMin = 0.0
betaMax = 10.0

beta_i = 0.5




lambda_s = np.exp(-beta_i*0.5) / ( 2.0*beta_i*np.sinh(beta_i*0.5) )
betaGrid = np.linspace(betaMin,betaMax,(betaMax-betaMin)/beta_i,endpoint=False)
continSAB = [0.0]*len(betaGrid)
for n,beta in enumerate(betaGrid):
    continSAB[n] = np.exp(-alpha*lambda_s) * 1.0/f(n) * (alpha*lambda_s)**n


plt.plot(betaGrid,continSAB,'ro')
plt.show()















