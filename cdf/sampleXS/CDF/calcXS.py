import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *
from rho import continRho,oscE,oscW
import random
import matplotlib.colors as colors
import matplotlib.cm as cmx
import seaborn as sns  # for nicer graphics

def prepPlot(alphas):
    plt.clf()
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    return scalarMap, plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)


def getFullSAB(alphas,betas,fullBetas,sabs): 
    fullSAB = [0.0]*(2*len(betas)-1)*len(alphas)
    for a in range(len(alphas)):
        for b in range(len(fullBetas)):
            bForLookup = abs(b-len(betas)+1)
            assert(abs(fullBetas[b]) == betas[bForLookup])
            fullSAB[a*len(fullBetas)+b] = sabs[a*len(betas)+bForLookup]
            if (fullBetas[b] > 0): fullSAB[a*len(fullBetas)+b] *= np.exp(-fullBetas[b])
    return fullSAB


def getAlphaMinMax(E,beta,kb,T,A):
    return ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T ),\
           ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )


def getAlphaMinMaxIndices(E,beta,kb,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    if aMin < alphas[0] or aMax > alphas[-1]:
        print("May want to have more alpha values")
    aMin_index, aMax_index = 0, 0
    for a in range(len(alphas)):
        if alphas[a] > aMin and aMin_index == 0: aMin_index = a
        if alphas[a] <= aMax:                    aMax_index = a
    return aMin_index,aMax_index


def getValidBetasRange(A,E,T,fullBetas):
    kb = 8.6173303e-5
    betaMin, betaMax = -E/(kb*T), 20.0
    if betaMin < fullBetas[0] or betaMax > fullBetas[-1]:
        print("May want to have more beta values")
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
            sabL = fullSAB[a*len(fullBetas)+b]
            sabR = fullSAB[(a+1)*len(fullBetas)+b]
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
    """
    plt.plot(fullBetas[bMin:bMax-1],eq14)
    plt.xlabel('beta'); plt.ylabel('beta PDF')
    plt.show()
    """
    return eq14



def calcBetaCDF(betas,eq14,bMin):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5 * \
                  (fullBetas[bMin+b]-fullBetas[bMin+b-1]) 
    return eq16





def getBetaCDF(A,E,T,fullSAB,fullBetas):
    bMin,bMax = getValidBetasRange(A,E,T,fullBetas)
    betaPDF = calcBetaPDF(A,E,T,fullSAB,fullBetas,bMin,bMax)
    betaCDF = calcBetaCDF(fullBetas,betaPDF,bMin)
    """
    plt.plot(fullBetas[bMin:bMax-1],betaCDF)
    plt.xlabel('beta'); plt.ylabel('beta CDF')
    plt.show()
    exit() 
    """
    return fullBetas[bMin:bMax-1],betaCDF


def binarySearch(vec,val):
    if len(vec) == 1: return 0
    if len(vec) == 2: return 1 if abs(vec[1]-val) < 1e-6 else 0
    if len(vec) == 3:
        return 0 if vec[0] <= val < vec[1] else \
               1 if vec[1] <= val < vec[2] else 2
    halfway = int(len(vec)*0.5)
    if abs(vec[halfway]-val) < 1e-6: return halfway
    return binarySearch(vec[:halfway],val) if val <= vec[halfway] else \
           binarySearch(vec[halfway:],val) + halfway



def sampleCDF(x,CDF):
    xsi = random.random()
    i = binarySearch(CDF,xsi)
    frac  = (xsi-CDF[i])/(CDF[i+1]-CDF[i])
    return i,frac*(x[i+1]-x[i])+x[i],xsi

def getAlphaCDF(nAlpha,nFullBetas,index,fullSAB,alphas,aMin):
    denominator = 0.0
    for a in range(nAlpha-1):
        sabL = fullSAB[a*nFullBetas+index]
        sabR = fullSAB[(a+1)*nFullBetas+index]
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    alphaPDF = [fullSAB[a*nFullBetas+index]/denominator \
                for a in range(nAlpha)]

    alphaCDF = [0.0]*nAlpha
    for a in range(1,nAlpha):
        alphaCDF[a] = alphaCDF[a-1]+(alphaPDF[a]+alphaPDF[a-1]) * \
                      0.5*(alphas[a]-alphas[a-1])
    return alphaCDF

A, T = 1.0, 296.0
kb = 8.6173303e-5

nAlpha, nBeta = 50, 100
alphas = list(np.linspace(0.001,3,nAlpha))
betas  = list(np.logspace(-6,2.1,nBeta))

nAlpha, nBeta = 50, 100
alphas = list(np.linspace(0.001,3,nAlpha))
betas  = list(np.logspace(-6,2.1,nBeta))

nAlpha, nBeta = 100, 1000
alphas = list(np.linspace(1.0,32.0,nAlpha))
alphas = list(np.linspace(0.0,32.0,nAlpha))
betas  = list(np.logspace(-6,1.3,nBeta))
print(kb*T)





useNJOY  = True
width    = None 
fullRedo = True
fullRedo = False  


SAB = getSAB(alphas,betas,T,continRho,useNJOY,fullRedo,width,oscE,oscW) 

"""
scalarMap, colorBar = prepPlot(alphas)
for a in range(0,nAlpha):
    plt.plot(betas,[SAB[a*len(betas)+b] for b in range(len(betas))],label=str("alpha = "+"%.3E"%alphas[a]),color=scalarMap.to_rgba(a),alpha=0.6)

plt.colorbar(colorBar).ax.set_ylabel('alpha')
plt.yscale('log')
plt.xlabel('beta'); plt.ylabel('S_n.sym(a,-b)')
plt.show()
"""

fullBetas = [-x for x in betas[1:]][::-1] + betas
fullSAB = getFullSAB(alphas,betas,fullBetas,SAB)

"""
scalarMap, colorBar = prepPlot(alphas)
for a in range(0,nAlpha):
    plt.plot(fullBetas,[fullSAB[a*len(fullBetas)+b] for b in range(len(fullBetas))],label=str("alpha = "+"%.3E"%alphas[a]),color=scalarMap.to_rgba(a),alpha=0.4)
plt.colorbar(colorBar).ax.set_ylabel('alpha')
plt.yscale('log')
plt.xlabel('beta'); plt.ylabel('S_n.sym(a,b)')
plt.show()
"""


N = int(5e6)

Energies = [0.0005,0.0253,0.2907,0.95,3.12]
#Energies = [3.12]
Energies = [0.0255]
xsiVals = []
betaVals = []

for i,E in enumerate(Energies):
    #break
    print("E = ",E)

    
    bMin,bMax = getValidBetasRange(A,E,T,fullBetas)
    CDF_beta_vals, betaCDF = getBetaCDF(A,E,T,fullSAB,fullBetas)

    E_out_vec = []

    xsi_Vals = []
    for n in range(N):
        if n%50000 == 0:
            print(n)
        index,beta,xsi = sampleCDF(CDF_beta_vals,betaCDF)
        xsiVals.append(xsi)
        betaVals.append(beta)

        """
        aMin, aMax = getAlphaMinMaxIndices(E,beta,kb,T,A,alphas)
        alphaCDF = getAlphaCDF(nAlpha,len(fullBetas),index,fullSAB,alphas,aMin)
        # If the valid alpha values for this chosen beta are not going to occur
        # we should just continue
        if (alphaCDF[aMax] < 1e-12): continue
        H = [None]*nAlpha
        for a in range(nAlpha):
            if   (alphas[a] <= alphas[aMin]): H[a] = 0
            elif (alphas[a] >= alphas[aMax]): H[a] = 1
            else: 
                H[a] = ( (alphaCDF[a]-alphaCDF[aMin]) / \
                         (alphaCDF[aMax]-alphaCDF[aMin]) )
        index,alpha = sampleCDF(alphas,H)
        """

        E_out_vec.append((beta*kb*T+E))
    plt.hist(E_out_vec,normed=True,bins=np.logspace(-4.5,0,80),alpha=0.3)
    print(min(E_out_vec),min(xsiVals),min(betaVals))
    #plt.hist(E_out_vec,normed=True,alpha=0.3)

plt.gca().set_xscale("log")
plt.yscale('log', nonposy='clip')
axes = plt.gca()
plt.xlabel("E' [eV]")
plt.ylabel("Probability")

axes.set_xlim([4e-5,1e0])
#axes.set_ylim([1e-3,1e3])


plt.show()
"""
"""












