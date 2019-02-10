import numpy as np
import matplotlib.pyplot as plt
from math import pi
from test09_sab import *

test9_alphas = alphas[:]

def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def calcAsym(alpha,beta):
    return (1.0/(4.0*pi*alpha)**0.5) * np.exp( - (alpha+beta)**2 / (4.0*alpha) )

def calcSym(alpha,beta):
    return np.exp(beta*0.5)*calcAsym(alpha,beta)

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax

n_beta = len(betas)
freeSAB = [0.0]*(2*len(betas)-1)*len(alphas)
fullBetas = [-x for x in betas[1:]][::-1] + betas
for a in range(len(alphas)):
    for b in range(len(fullBetas)):
        freeSAB[a*len(fullBetas)+b] = (calcAsym(alphas[a],fullBetas[b]))

fullSAB = [0.0]*(2*len(betas)-1)*len(alphas)
for a in range(len(alphas)):
    for b in range(len(fullBetas)):
        b_for_positive_betas = abs(b-len(betas)+1)
        beta = betas[b_for_positive_betas]
        assert(abs(fullBetas[b]) == beta)
        # This creates the symmetric SAB
        fullSAB[a*len(fullBetas)+b] = getSABval(sab,a,b_for_positive_betas,len(betas))
        #fullSAB[a*len(fullBetas)+b] = (calcSym(alphas[a],fullBetas[b]))
        # We make it asymmetric to make sure it matches eq14
        fullSAB[a*len(fullBetas)+b] *= np.exp(-fullBetas[b]*0.5)



def getValidBetasRange(A,E,T,fullBetas):
    kb = 8.6173303e-5
    betaMin = -E/(kb*T)
    betaMax = 20.0

    bMin = bMax  = len(fullBetas)
    for b in range(len(fullBetas)):
        if bMin == len(fullBetas) and betaMin <= fullBetas[b]:
                bMin = b
        if fullBetas[b] >= betaMax:
            bMax = b
            break
    return bMin,bMax



def calcIntegralAcrossAlpha(A,E,T,fullBetas,b):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,fullBetas[b],kb,T,A)
    denominator = 0.0
    for a in range(len(alphas)-1):
        if aMin <= alphas[a] <= aMax:
            #sabL = getSABval(freeSAB,a,b,len(fullBetas))
            #sabR = getSABval(freeSAB,a+1,b,len(fullBetas))
            sabL = getSABval(fullSAB,a,b,len(fullBetas))
            sabR = getSABval(fullSAB,a+1,b,len(fullBetas))

            denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    return denominator

def calcBetaPDF(A,E,T,fullBetas,bMin,bMax):

    denominator = 0.0
    delta0 = fullBetas[bMin+1]-fullBetas[bMin]
    deltaN = fullBetas[bMax-1]-fullBetas[bMax-2]
    denominator = calcIntegralAcrossAlpha(A,E,T,fullBetas,bMin)*delta0*0.5 +\
                  calcIntegralAcrossAlpha(A,E,T,fullBetas,bMax-1)*deltaN*0.5
    for b in range(bMin+1,bMax-1):
        deltaL = fullBetas[b]-fullBetas[b-1]
        deltaR = fullBetas[b+1]-fullBetas[b]
        denominator += calcIntegralAcrossAlpha(A,E,T,fullBetas,b) * (deltaL+deltaR) * 0.5

    eq14 = [ calcIntegralAcrossAlpha(A,E,T,fullBetas,b) / denominator \
           for b in range(bMin,bMax-1)]
    return eq14



def calcBetaCDF(betas,eq14):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5*(fullBetas[b]-fullBetas[b-1]) 
    return eq16


def PDF_CDF_at_various_temperatures(A,E,temps,fullBetas):
    betas_vec = []
    eq14_vec = []
    eq16_vec = []
    for T in temps:
        bMin,bMax = getValidBetasRange(A,E,T,fullBetas)
        eq14 = calcBetaPDF(A,E,T,fullBetas,bMin,bMax)
        eq16 = calcBetaCDF(fullBetas,eq14)
        betas_vec.append(fullBetas[bMin:bMax-1])
        eq14_vec.append(eq14)
        eq16_vec.append(eq16)
    for i in range(len(temps)):
        plt.plot(betas_vec[i],eq14_vec[i])
    plt.show()
    for i in range(len(temps)):
        plt.plot(betas_vec[i],eq16_vec[i])
    plt.show()


A = 18.0
E = 1.0 
temps = [300.0,475.0,650.0,825.0,1000.0]   # Water
temps = [300.0]   # Water
PDF_CDF_at_various_temperatures(A,E,temps,fullBetas)

