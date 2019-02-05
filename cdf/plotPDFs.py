import matplotlib.pyplot as plt
from plotHelp import *
from test09_sab import *
import numpy as np


def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]


fullBetas = [-x for x in betas[::-1][:-1]]+betas
fullSAB = [[0.0]*len(fullBetas) for a in range(len(alphas))]
for a in range(n_alpha):
    for b in range(-n_beta+1,n_beta):
        fullSAB[a][n_beta+b-1] = getSABval(sab,a,abs(b),n_beta)


plt.plot(fullBetas,fullSAB[0])

# find denominator
eq14_denom = [0.0]*len(fullBetas)
for b in range(len(fullBetas)-1):
    exp_term = np.exp(-fullBetas[b]*0.5)
    val_b_left  = 0.0
    val_b_right = 0.0
    for a in range(n_alpha-1):
        val_b_left  += (fullSAB[a][b]+fullSAB[a+1][b])  *0.5*(alphas[a+1]-alphas[a])
        val_b_right += (fullSAB[a][b]+fullSAB[a+1][b+1])*0.5*(alphas[a+1]-alphas[a])
    eq14_denom[b] += exp_term*(val_b_left+val_b_right)*0.5**(fullBetas[b+1]-fullBetas[b])


eq14 = [0.0]*len(fullBetas)
for b in range(len(fullBetas)-1):
    exp_term = np.exp(-fullBetas[b]*0.5)
    val = 0.0
    for a in range(n_alpha-1):
        val += (fullSAB[a][b]+fullSAB[a+1][b]) * 0.5 * (alphas[a+1]-alphas[a])
    eq14[b] = exp_term * val / eq14_denom[b]




plt.plot(fullBetas,eq14_denom)
plt.plot(fullBetas,eq14_denom,'o')
plt.plot(fullBetas,eq14)
plt.plot(fullBetas,eq14,'o')
plt.yscale('log')

plt.show()



