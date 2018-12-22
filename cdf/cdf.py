import matplotlib.pyplot as plt
#from results_test_09_300K import *
#from results_test_09_300K_with_negative_beta import *
#from results_symBeta import *
from resultsFreeGas import *
#from resultsMacFarlane import *
#from resultsMacFarlane2016_296K import *
from math import exp
from math import pi
import numpy as np

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax


def calcPDF(betaVec,E,kb,T,A,alphas,sab,nbeta,sabIsSym,sign):
    eq14 = []
    for b,beta in enumerate(betaVec):
        aMin, aMax = getAlphaMinMax(E,sign*beta,kb,T,A)
        sabVals = [sab[a*nbeta+b] for a,alpha in enumerate(alphas) if (aMin<=alpha<=aMax)]
        validAlphas = [ a for a in alphas if (a >= aMin and a <= aMax) ]
        eq14 = eq14 + [exp(-0.5*sign*beta)*np.trapz(sabVals,x=validAlphas) \
                          if sabIsSym else np.trapz(sabVals,x=validAlphas)]
    return eq14



A = 1.0  
E = 1    # 1 eVc
T = 296  # 300 K
kb = 8.6173303e-5

bMin = -E/(kb*T)
bMax = 20.0

sabIsSym = True
sabIsSym = False
#sab = sab_water
nbeta = len(betas)

#print(alphas[0],betas[0])
#print(sab[0*nbeta+0])


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )

def calcSym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return exp(beta*0.5)*calcAsym(alpha,beta)

#print(calcSym(alphas[0],betas[0]))
#print(calcAsym(alphas[0],betas[0]))


# Define valid beta vectors for +Beta, -Beta
betaPosVec = [ b for b in betas if b <= bMax      ]
betaNegVec = [ b for b in betas if b <= abs(bMin) ]

# Calculate respective PDFs
eq14Pos = calcPDF(betaPosVec,E,kb,T,A,alphas,sab,nbeta,sabIsSym,+1)
eq14Neg = calcPDF(betaNegVec,E,kb,T,A,alphas,sab,nbeta,sabIsSym,-1)

# Combine the positive and negative beta regions and plot
eq14 = eq14Neg[::-1] + eq14Pos
betaTotal = [-b for b in betaNegVec[::-1]] + betaPosVec
invIntegral = 1.0/np.trapz(eq14, x=betaTotal)
eq14 = [x * invIntegral for x in eq14]
#print(eq14)

f = plt.figure(1); plt.plot(betaTotal,eq14); f.show(); input()


