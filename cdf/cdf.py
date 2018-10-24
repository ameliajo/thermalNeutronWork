import matplotlib.pyplot as plt
#from results_test_09_300K import *
#from results_test_09_300K_with_negative_beta import *
#from results_symBeta import *
from resultsFreeGas import *
from math import exp
import numpy as np

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax


def getEq14(beta,E,T,Sab,A,alphaVec,lenBeta,b):
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    validAlphaIndexes = [ i for i in range(len(alphaVec)) if \
                   (alphaVec[i] >= aMin and alphaVec[i] <= aMax) ]


    sabVec = [Sab[i*lenBeta+b] for i in validAlphaIndexes]
    validAlpha = alphaVec[validAlphaIndexes[0]:validAlphaIndexes[-1]+1]

    numerator = exp(-beta/2) * np.trapz([Sab[i*lenBeta+b] for i in validAlphaIndexes],x=alphaVec[validAlphaIndexes[0]:validAlphaIndexes[-1]+1])

    return numerator 



A = 1.0  
E = 1    # 1 eVc
T = 1000  # 300 K
kb = 8.6173303e-5

bMin = -E/(kb*T)
bMax = 20.0


sabIsSym = True
sab = sab_water

eq14Pos = []
posBMin = 0.0
posBMax = bMax

betaPosVec = [ b for b in betas if (b >= posBMin and b <= posBMax) ]
for b,beta in enumerate(betaPosVec):
    bIndex = betas.index(betaPosVec[b])
    aMin, aMax = getAlphaMinMax(E,beta,kb,T,A)
    validAlphas = [ a for a in alphas if (a >= aMin and a <= aMax) ]
    
    sabVals = [ sab[alphas.index(alpha)*len(betas)+bIndex] for alpha in validAlphas ]
    if sabIsSym: eq14Pos.append(exp(-0.5*beta)*np.trapz(sabVals,x=validAlphas))
    else       : eq14Pos.append(               np.trapz(sabVals,x=validAlphas))
 
eq14Neg = []
negBMin = 0.0
negBMax = abs(bMin)

betaNegVec = [ abs(b) for b in betas if (b >= negBMin  and b <= negBMax) ]
for b,beta in enumerate(betaNegVec):
    bIndex = betas.index(betaNegVec[b])
    aMin, aMax = getAlphaMinMax(E,-beta,kb,T,A)
    validAlphas = [ a for a in alphas if (a >= aMin and a <= aMax) ]
    sabVals = [ sab[alphas.index(alpha)*len(betas)+bIndex] for alpha in validAlphas ]

    if sabIsSym: eq14Neg.append(exp(0.5*beta)*np.trapz(sabVals,x=validAlphas))
    else       : eq14Neg.append(              np.trapz(sabVals,x=validAlphas))


       
eq14 = eq14Neg[::-1] + eq14Pos
betaTotal = [-b for b in betaNegVec[::-1]] + betaPosVec
invIntegral = 1.0/np.trapz(eq14, x=betaTotal)
eq14 = [x * invIntegral for x in eq14]
print(eq14)

f = plt.figure(1); plt.plot(betaTotal,eq14); f.show(); input()


