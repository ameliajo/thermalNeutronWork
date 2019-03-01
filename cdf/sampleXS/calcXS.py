import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *
from rho import continRho,oscE,oscW
from random import random

def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]


def getFullSAB(alphas,betas,fullBetas,sabs): 
    fullSAB = [0.0]*(2*len(betas)-1)*len(alphas)
    for a in range(len(alphas)):
        for b in range(len(fullBetas)):
            bForLookup = abs(b-len(betas)+1)
            assert(abs(fullBetas[b]) == betas[bForLookup])
            # This creates the symmetric SAB
            fullSAB[a*len(fullBetas)+b] = getSABval(sabs,a,bForLookup,len(betas))
            # We make it asymmetric to make sure it matches eq14
            if (fullBetas[b] > 0): fullSAB[a*len(fullBetas)+b] *= np.exp(-fullBetas[b])
    return fullSAB


def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax


def getAlphaMinMaxIndices(E,beta,kb,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    aMin_index, aMax_index = 0, 0
    for a in range(len(alphas)):
        if alphas[a] > aMin and aMin_index == 0: aMin_index = a
        if alphas[a] <= aMax:                    aMax_index = a
    return aMin_index,aMax_index





def getValidBetasRange(A,E,T,fullBetas):
    kb = 8.6173303e-5
    betaMin, betaMax = -E/(kb*T), 20.0
    bMin = bMax  = len(fullBetas)
    for b in range(len(fullBetas)):
        if bMin == len(fullBetas) and betaMin <= fullBetas[b]:
            bMin = b
        if fullBetas[b] >= betaMax:
            return bMin,b
    return bMin,bMax


def intSABda(A,E,T,fullSAB,fullBetas,b):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,fullBetas[b],kb,T,A)
    denominator = 0.0
    for a in range(len(alphas)-1):
        if aMin <= alphas[a] <= aMax:
            sabL = getSABval(fullSAB,a,b,len(fullBetas))
            sabR = getSABval(fullSAB,a+1,b,len(fullBetas))
            denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    return denominator

def calcBetaPDF(A,E,T,fullSAB,fullBetas,bMin,bMax):
    delta0 = fullBetas[bMin+1]-fullBetas[bMin]
    deltaN = fullBetas[bMax-1]-fullBetas[bMax-2]
    denominator = intSABda(A,E,T,fullSAB,fullBetas,bMin)*delta0*0.5 +\
                  intSABda(A,E,T,fullSAB,fullBetas,bMax-1)*deltaN*0.5
    for b in range(bMin+1,bMax-1):
        deltaL = fullBetas[b]-fullBetas[b-1]
        deltaR = fullBetas[b+1]-fullBetas[b]
        denominator += intSABda(A,E,T,fullSAB,fullBetas,b) * (deltaL+deltaR) * 0.5

    eq14 = [ intSABda(A,E,T,fullSAB,fullBetas,b) / denominator \
           for b in range(bMin,bMax-1)]
    return eq14



def calcBetaCDF(betas,eq14,bMin):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5 * \
                  (fullBetas[bMin+b]-fullBetas[bMin+b-1]) 
    return eq16





def PDF_CDF_at_various_temperatures(A,E,T,fullSAB,fullBetas):
    bMin,bMax = getValidBetasRange(A,E,T,fullBetas)
    eq14 = calcBetaPDF(A,E,T,fullSAB,fullBetas,bMin,bMax)
    eq16 = calcBetaCDF(fullBetas,eq14,bMin)
    return fullBetas[bMin:bMax-1],eq16




E, A, T = 1.0, 1.0, 296.0
kb = 8.6173303e-5

n_alpha, n_beta = 500, 100
alphas = list(np.linspace(0.01,60,n_alpha))
betas  = list(np.linspace(0.00,40,n_beta))

useNJOY  = True
fullRedo = False
width    = None 


fullBetas = [-x for x in betas[1:]][::-1] + betas
SAB = getSAB(alphas,betas,T,continRho,useNJOY,fullRedo,width,oscE,oscW) 
fullSAB = getFullSAB(alphas,betas,fullBetas,SAB)

relevantBetas, betaCDF = PDF_CDF_at_various_temperatures(A,E,T,fullSAB,fullBetas)

N = 1000
for n in range(N):
    xsi_1 = random()
    index = 0
    beta = 0.0
    for i in range(len(relevantBetas)-1):
        if betaCDF[i] <= xsi_1 < betaCDF[i+1]:
            index = i
            frac  = (xsi_1-betaCDF[i])/(betaCDF[i+1]-betaCDF[i])
            beta = frac*(betaCDF[i+1]-betaCDF[i])+relevantBetas[i]
            break


    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(fullSAB,a,index,len(fullBetas))
        sabR = getSABval(fullSAB,a+1,index,len(fullBetas))
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    alphaPDF = [getSABval(fullSAB,a,index,len(fullBetas))/denominator \
                for a in range(len(alphas))]
    
    alphaCDF = [0.0]*len(alphaPDF)
    for a in range(1,len(alphaPDF)):
        alphaCDF[a] = alphaCDF[a-1]+(alphaPDF[a]+alphaPDF[a-1]) * \
                      0.5*(alphas[a]-alphas[a-1])

    aMin_index, aMax_index = getAlphaMinMaxIndices(E,beta,kb,T,A,alphas)
    alphasToPlot, H = [], []
    for a in range(aMin_index,aMax_index+1):
        alphasToPlot.append(alphas[a])
        H.append( (alphaCDF[a]-alphaCDF[aMin_index]) / \
                  (alphaCDF[aMax_index]-alphaCDF[aMin_index]) )
    

    xsi_2 = random()
    alpha = 0.0
    for i in range(len(alphasToPlot)-1):
        if alphas[i] >= alphasToPlot[0] and alphas[i] <= alphasToPlot[-1]:
            if H[i] <= xsi_2 < H[i+1]:
                frac  = (xsi_2-H[i])/(H[i+1]-H[i])
                alpha = frac*(H[i+1]-H[i])+alphasToPlot[i]
                break



    plt.plot(alpha,beta,'ro')











#plt.plot(relevantBetas,betaCDF)
plt.xlabel('alpha')
plt.ylabel('beta')
plt.show()











