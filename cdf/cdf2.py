import matplotlib.pyplot as plt
from sabFreeGasSym   import *
from sabFreeGas      import *
from alphaFreeGas    import *
from betaFreeGas     import *
import numpy as np
from math import pi
from math import exp

sab = sabSym[:]

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )


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
    return exp(beta*0.5)*calcAsym(alpha,beta)


A = 1.0 
E = 1      # 1 eVc
T = 1000    # 300 K
kb = 8.6173303e-5

bMin, bMax = -E/(kb*T), 20.0
sabIsSym = True



eq14Pos = []
posBMin = 0.0
posBMax = bMax

betaPosVec = [ b for b in betas if (b >= posBMin and b <= posBMax) ]
for b,beta in enumerate(betaPosVec):
    bIndex = betas.index(betaPosVec[b])
    aMin, aMax = getAlphaMinMax(E,beta,kb,T,A)
    validAlphas = [ a for a in alphas[bIndex] if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas[bIndex].index(alpha)
        if     sabIsSym and (sab[bIndex][aIndex]-calcSym(alpha,beta))  > 1e-15: print("OH NO")
        if not sabIsSym and (sab[bIndex][aIndex]-calcAsym(alpha,beta)) > 1e-15: print("OH NO")
        sabVals.append(sab[bIndex][aIndex])
    if sabIsSym: eq14Pos.append(exp(-0.5*beta)*np.trapz(sabVals,x=validAlphas))
    else       : eq14Pos.append(               np.trapz(sabVals,x=validAlphas))
    print(beta,eq14Pos[-1])
 
eq14Neg = []
negBMin = 0.0
negBMax = abs(bMin)

betaNegVec = [ abs(b) for b in betas if (b >= negBMin  and b <= negBMax) ]
for b,beta in enumerate(betaNegVec):
    bIndex = betas.index(-betaNegVec[b])
    aMin, aMax = getAlphaMinMax(E,-beta,kb,T,A)
    validAlphas = [ a for a in alphas[bIndex] if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas[bIndex].index(alpha)
        if     sabIsSym and (sab[bIndex][aIndex]-calcSym(alpha,-beta))  > 1e-15: print("OH NO")
        if not sabIsSym and (sab[bIndex][aIndex]-calcAsym(alpha,-beta)) > 1e-15: print("OH NO")
        sabVals.append(sab[bIndex][aIndex])
    if sabIsSym: eq14Neg.append(exp(0.5*beta)*np.trapz(sabVals,x=validAlphas))
    else       : eq14Neg.append(              np.trapz(sabVals,x=validAlphas))


       
eq14Neg.reverse()
betaNegVec.reverse()
eq14 = eq14Neg + eq14Pos
betaTotal = [-b for b in betaNegVec] + betaPosVec
invIntegral = 1.0/np.trapz(eq14, x=betaTotal)
eq14 = [x * invIntegral for x in eq14]
#print(eq14)
f = plt.figure(1)
plt.plot(betaTotal,eq14)
f.show()
#hold = input()


if __name__ == "__main__":
    import doctest
    doctest.testmod()
        
