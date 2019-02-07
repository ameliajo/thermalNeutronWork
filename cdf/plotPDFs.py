import matplotlib.pyplot as plt
from plotHelp import *
#from test09_sab import *
import numpy as np
from math import pi

def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return (1.0/(4.0*pi*alpha)**0.5) * np.exp( - (alpha+beta)**2 / (4.0*alpha) )


def calcSym(alpha,beta):
    """
    >>> abs(calcSym(1.0,2.0)-calcSym(1.0,-2.0)) < 1e-15
    True
    >>> abs(calcSym(5.0,12.0)-calcSym(5.0,-12.0)) < 1e-15
    True
    >>> abs(calcSym(3.5,5.0)-calcSym(2.5,-5.0)) < 1e-15
    False
    >>> abs(calcSym(2.5,4.0)-calcSym(2.5,-5.0)) < 1e-15
    False
    """
    if alpha == 0: alpha = 0.0001
    return np.exp(beta*0.5)*calcAsym(alpha,beta)

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax



A = 18.0 
E = 1.0      # 1 eV
T = 300    # 300 K
kb = 8.6173303e-5

bMin, bMax = -E/(kb*T), 20.0

alphas = list(np.linspace(0.1,25,500))
fullBetas = list(np.linspace(bMin,bMax,100))
n_alpha = len(alphas)
n_beta = len(fullBetas)


fullSAB = [[0.0]*len(fullBetas) for a in range(len(alphas))]
for a in range(n_alpha):
    for b in range(n_beta):
        fullSAB[a][b] = calcSym(a,b)


# find denominator
eq14_denom = [0.0]*len(fullBetas)
for b in range(len(fullBetas)-1):
    exp_term = np.exp(-fullBetas[b]*0.5)
    val_b_left, val_b_right = 0.0, 0.0
    a_min, a_max = getAlphaMinMax(1.0,abs(fullBetas[b]),8.6173303e-5,1000.0,18)
    for a in range(n_alpha-1):
        if ( a_min < alphas[a] and alphas[a+1] < a_max):
            val_b_left  += (fullSAB[a][b]+fullSAB[a+1][b])  *0.5*(alphas[a+1]-alphas[a])
            val_b_right += (fullSAB[a][b]+fullSAB[a+1][b+1])*0.5*(alphas[a+1]-alphas[a])
    eq14_denom[b] += exp_term*(val_b_left+val_b_right)*0.5**(fullBetas[b+1]-fullBetas[b])


eq14 = [0.0]*len(fullBetas)
for b in range(len(fullBetas)-1):
    exp_term = np.exp(-fullBetas[b]*0.5)
    val = 0.0
    a_min, a_max = getAlphaMinMax(1.0,abs(fullBetas[b]),8.6173303e-5,1000.0,18)
    for a in range(n_alpha-1):
        if ( a_min < alphas[a] and alphas[a+1] < a_max):
            val += (fullSAB[a][b]+fullSAB[a+1][b]) * 0.5 * (alphas[a+1]-alphas[a])
    eq14[b] = exp_term * val / eq14_denom[b]




plt.plot(fullBetas[:-1],eq14[:-1])
#plt.plot(fullBetas[:-1],eq14[:-1],'o')
#plt.yscale('log')
plt.show()



