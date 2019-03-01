import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *
from rho import continRho,oscE,oscW
import random

import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20c')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('energies')
    #plt.title(title)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    #ax.set_facecolor('xkcd:very light blue') # off white
    #plt.yscale('log')
    plt.show()





def getSABval(sab,a,b,nBeta):
    return sab[a*nBeta+b]


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


def sampleCDF(x,CDF):
    xsi = random.random()
    for i in range(len(x)-1):
        if CDF[i] <= xsi < CDF[i+1]:
            frac  = (xsi-CDF[i])/(CDF[i+1]-CDF[i])
            val = frac*(CDF[i+1]-CDF[i])+x[i]
            return i,val



def getAlphaCDF(nAlpha,nFullBetas,index,fullSAB,alphas):
    denominator = 0.0
    for a in range(nAlpha-1):
        sabL = getSABval(fullSAB,a,index,nFullBetas)
        sabR = getSABval(fullSAB,a+1,index,nFullBetas)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    alphaPDF = [getSABval(fullSAB,a,index,nFullBetas)/denominator \
                for a in range(nAlpha)]
    alphaCDF = [0.0]*nAlpha
    for a in range(1,nAlpha):
        alphaCDF[a] = alphaCDF[a-1]+(alphaPDF[a]+alphaPDF[a-1]) * \
                      0.5*(alphas[a]-alphas[a-1])
    return alphaCDF




A, T = 1.0, 296.0
kb = 8.6173303e-5

nAlpha, nBeta = 500, 100
alphas = list(np.linspace(0.01,60,nAlpha))
betas  = list(np.linspace(0.00,40,nBeta))

useNJOY  = True
fullRedo = False
width    = None 


fullBetas = [-x for x in betas[1:]][::-1] + betas
SAB = getSAB(alphas,betas,T,continRho,useNJOY,fullRedo,width,oscE,oscW) 
fullSAB = getFullSAB(alphas,betas,fullBetas,SAB)

N = 500


Energies = [0.01,0.1,0.5,1.0,5.0,6.0,8.0,10.0,12.0]
scalarMap, colorBar = prepPlot(Energies)
for i,E in enumerate(Energies):
    print(E)

    relevantBetas, betaCDF = PDF_CDF_at_various_temperatures(A,E,T,fullSAB,fullBetas)

    for n in range(N):
    
        index,beta = sampleCDF(relevantBetas,betaCDF)
    
        alphaCDF = getAlphaCDF(nAlpha,len(fullBetas),index,fullSAB,alphas)
    
        aMin, aMax = getAlphaMinMaxIndices(E,beta,kb,T,A,alphas)

        H = [None]*nAlpha
        for a in range(nAlpha):
            if   (alphas[a] < alphas[aMin]): H[a] = 0
            elif (alphas[a] > alphas[aMax]): H[a] = 1
            else:
                H[a] = ( (alphaCDF[a]-alphaCDF[aMin]) / \
                         (alphaCDF[aMax]-alphaCDF[aMin]) )
    
        index,alpha = sampleCDF(alphas,H)

        plt.plot(alpha,beta,marker='o',markersize=1,c=scalarMap.to_rgba(i))





#plt.plot(relevantBetas,betaCDF)
plt.xlabel('alpha')
plt.ylabel('beta')
finishPlotting(colorBar,"")











