import matplotlib.pyplot as plt
from sabFreeGas   import *
from alphaFreeGas import *
from betaFreeGas  import *
import numpy as np
from math import pi
from math import exp


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )



for b in range(len(betas)):
    beta = betas[b]
    alphaVec = alphas[b] 
    for a in range(len(alphaVec)):
        alpha = alphaVec[a]
        if abs(sab[b][a]-calcAsym(alpha,beta))>1e-15: 
            print(beta,alpha,sab[b][a]-calcAsym(alpha,beta))



def getEq14(beta,E,T,Sab,A,alphas,lenBeta,b):
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aIndices = [ i for i in range(len(alphas)) if (alphas[i] >= aMin and alphas[i] <= aMax) ]

    sabVec = [sab[b][a] for a in aIndices]
    validAlpha = alphaVec[aIndices[0]:aIndices[-1]+1]
    return np.trapz(sabVec,x=validAlpha)

    #numerator = exp(-beta/2) * np.trapz([sab[b][a] for i in aIndices],x=alphaVec[aIndices[0]:aIndices[-1]+1])






A = 1.0 
E = 1      # 1 eVc
T = 1000    # 300 K
kb = 8.6173303e-5

bMin, bMax = -E/(kb*T), 20.0


"""

eq14, eq14Neg = [],[]
eq14Neg = []
for b,beta in enumerate(betas):

    alphaVec = alphas[b] 
    eq14.append(getEq14(beta,E,T,sab,A,alphaVec,len(betas),b))
    eq14Neg.append(getEq14(-beta,E,T,sab,A,alphaVec,len(betas),b))
    print(eq14Neg[0])
    break 
    #eq14.append(getEq14(beta,E,T,sab_water,A,alphas,len(betas),b))
    #eq14Neg.append(getEq14(-beta,E,T,sab_water,A,alphas,len(betas),b))

#eq14Neg.reverse()
#eq14Full = eq14Neg + eq14
#betaFull = 

"""


eq14Pos = []
posBMin = 0.0
posBMax = bMax

betaPosVec = [ b for b in betas if (b >= posBMin and b <= posBMax) ]
for b in range(len(betaPosVec)):
    beta = betaPosVec[b]
    bIndex = betas.index(betaPosVec[b])
    alphaVec = alphas[bIndex] 
    aMin = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    validAlphas = [ a for a in alphaVec if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas[bIndex].index(alpha)
        if (sab[bIndex][aIndex]-calcAsym(alpha,betaPosVec[b])) > 1e-15: print("OH NO")
        sabVals.append(sab[bIndex][aIndex])
    numerator = np.trapz(sabVals,x=validAlphas)
    eq14Pos.append(numerator)

#print(eq14Pos)

 
eq14Neg = []
negBMin = 0.0
negBMax = abs(bMin)

betaNegVec = [ abs(b) for b in betas if (b >= negBMin  and b <= negBMax) ]
for b in range(len(betaNegVec)):
    beta = betaNegVec[b]
    bIndex = betas.index(-betaNegVec[b])
    alphaVec = alphas[bIndex] 
    aMin = ( (E)**0.5 - (E+(-beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+(-beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    validAlphas = [ a for a in alphaVec if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas[bIndex].index(alpha)
        if (sab[bIndex][aIndex]-calcAsym(alpha,-betaNegVec[b])) > 1e-15: print("OH NO")
        sabVals.append(sab[bIndex][aIndex])
    numerator = np.trapz(sabVals,x=validAlphas)
    eq14Neg.append(numerator)

       
eq14Neg.reverse()
betaNegVec.reverse()
print(eq14Neg)
eq14 = eq14Neg + eq14Pos
betaTotal = [-b for b in betaNegVec] + betaPosVec
print(betaTotal)
invIntegral = 1.0/np.trapz(eq14, x=betaTotal)
eq14 = [x * invIntegral for x in eq14]
f = plt.figure(1)
plt.plot(betaTotal,eq14)
f.show()
hold = input()


        
