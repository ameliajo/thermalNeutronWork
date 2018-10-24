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

print(bMin)

sabIsSym = True
sab = sab_water

eq14Pos = []
posBMin = 0.0
posBMax = bMax

betaPosVec = [ b for b in betas if (b >= posBMin and b <= posBMax) ]
print(betaPosVec)
for b,beta in enumerate(betaPosVec):
    bIndex = betas.index(betaPosVec[b])
    aMin, aMax = getAlphaMinMax(E,beta,kb,T,A)
    validAlphas = [ a for a in alphas if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas.index(alpha)
        sabVals.append(sab[aIndex*len(betas)+bIndex])
    if sabIsSym: eq14Pos.append(exp(-0.5*beta)*np.trapz(sabVals,x=validAlphas))
    else       : eq14Pos.append(               np.trapz(sabVals,x=validAlphas))
    #print(beta,eq14Pos[-1])
 
eq14Neg = []
negBMin = 0.0
negBMax = abs(bMin)

betaNegVec = [ abs(b) for b in betas if (b >= negBMin  and b <= negBMax) ]
for b,beta in enumerate(betaNegVec):
    bIndex = betas.index(betaNegVec[b])
    aMin, aMax = getAlphaMinMax(E,-beta,kb,T,A)
    validAlphas = [ a for a in alphas if (a >= aMin and a <= aMax) ]
    
    sabVals = []
    for alpha in validAlphas:
        aIndex = alphas.index(alpha)
        sabVals.append(sab[aIndex*len(betas)+bIndex])
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
hold = input()


"""
eq14 = []
eq14Neg = []
for b,beta in enumerate(betas):
    eq14.append(getEq14(beta,E,T,sab_water,A,alphas,len(betas),b))
    eq14Neg.append(getEq14(-beta,E,T,sab_water,A,alphas,len(betas),b))


#invTotal = 1.0/(sum(eq14)+sum(eq14Neg))
b2 = betas[:]; b2.reverse()
fullbetas = [-x for x in b2] + betas

invTotal = 1.0/np.trapz(eq14Neg+eq14,x=fullbetas)
fullbetas = [x * invTotal for x in fullbetas]
totaleq14 = eq14Neg+eq14
plt.plot(fullbetas,totaleq14,'ro-',markersize=2)
print(fullbetas)
f=plt.figure(1)

"""
#eq14    = [ eq14Val * invTotal for eq14Val in eq14    ]
#eq14Neg = [ eq14Val * invTotal for eq14Val in eq14Neg ]
#plt.plot(betas,eq14,'ro-',markersize=2)
#plt.plot([-b for b in betas],eq14Neg,'ro-',markersize=2)
#plt.show()
#f.show()
#input()


