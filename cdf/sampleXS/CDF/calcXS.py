import numpy as np
from math import pi
#from getSAB import *
#from rho import continRho,oscE,oscW
import random
from plotHelp import *
from sab import SAB


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
    if betaMin < fullBetas[0] or betaMax > 1.1*fullBetas[-1]:
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
    integral_SAB_da = [intSABda(A,E,T,fullSAB,fullBetas,b) for b in range(bMin,bMax-1)]
    denominator = np.trapz(integral_SAB_da,x=fullBetas[bMin:bMax-1])
    # This is Eq.14 from Pavlou and Ji paper:
    #     On-the-fly sampling of temperature-dependent thermal neutron 
    #     scattering data for Monte Carlo simulations
    return [ val / denominator for val in integral_SAB_da ]



def calcBetaCDF(betas,eq14,bMin):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5 * \
                  (fullBetas[bMin+b]-fullBetas[bMin+b-1]) 
    return eq16




def getBetaCDF(A,E,T,fullSAB,fullBetas):
    # These are the indices of the min/max allowed beta values
    bMin,bMax = getValidBetasRange(A,E,T,fullBetas)
    betaPDF = calcBetaPDF(A,E,T,fullSAB,fullBetas,bMin,bMax)
    betaCDF = calcBetaCDF(fullBetas,betaPDF,bMin)
    #plotBetaCDF(fullBetas,betaCDF,bMin,bMax)
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


def findLocation(vec,val):
    for i in range(len(vec)-1):
        if vec[i] <= val and val < vec[i+1]:
            return i


def binarySearch2(vec,val):
    if len(vec) == 1:
        return 0
    halfway = int(len(vec)*0.5)
    if val < vec[halfway]:
        return binarySearch2(vec[:halfway],val)
    else:
        return halfway+binarySearch2(vec[halfway:],val)






def sampleCDF(x,CDF):
    xsi = random.random()
    i = binarySearch2(CDF,xsi)
    assert(CDF[i]<=xsi<=CDF[i+1])
    frac  = (xsi-CDF[i])/(CDF[i+1]-CDF[i])
    return i,frac*(x[i+1]-x[i])+x[i]

def getAlphaCDF(nAlpha,nFullBetas,index,fullSAB,alphas):

    alphaPDF = [fullSAB[a*nFullBetas+index] for a in range(nAlpha)]
    denom = np.trapz(alphaPDF,x=alphas)
    alphaPDF = [x/denom for x in alphaPDF]

    alphaCDF = [0.0]*nAlpha
    for a in range(1,nAlpha):
        alphaCDF[a] = alphaCDF[a-1]+(alphaPDF[a]+alphaPDF[a-1]) * \
                      0.5*(alphas[a]-alphas[a-1])
    #plotAlphaCDF(alphas,alphaCDF)
    return alphaCDF

A, T = 1.0, 296.0
kb = 8.6173303e-5

nAlpha, nBeta = 100, 1000
alphas = list(np.linspace(0.0,32.0,nAlpha))
betas  = list(np.logspace(-6,1.3,nBeta))

useNJOY  = True
width    = None 
#fullRedo = True
fullRedo = False  


## This reads in the S(a,b) from NJOY. The S(a,b) that LEAPR returns is the 
## non-symmetric form, and they give the -b side of it. 
## S_sym(a,b)   = exp(b/2) * S_n.sym(a,b)
## S_n.sym(a,b) = exp(-b)  * S_n.sym(a,-b)
# SAB = getSAB(alphas,betas,T,continRho,useNJOY,fullRedo,width,oscE,oscW) 

#with open('sab.py', 'w') as f:
#    f.write("sab = [")
#    for sabVal in SAB[:-1]:
#        f.write("%s, " % sabVal)
#    f.write("%s ]" % SAB[-1])
plotSAB(alphas,betas,SAB,'S_n.sym(a,-b)')
exit()

## Reflect beta so that we have all positive and negative betas
## Fill out full S(a,b_ so that it's defined for positive and negative betas
fullBetas = [-x for x in betas[1:]][::-1] + betas
fullSAB = getFullSAB(alphas,betas,fullBetas,SAB)
# plotSAB(alphas,fullBetas,fullSAB,'S_n.sym(a,b)')


# How many neutrons do you want to sample for each energy?
N = int(5e5)

# Initial neutron energy
Energies = [0.0005,0.0253,0.2907,0.95,3.12]

for i,E in enumerate(Energies):
    # Get beta CDF
    CDF_beta_vals, betaCDF = getBetaCDF(A,E,T,fullSAB,fullBetas)
    E_out_vec = []

    for n in range(N):
        # Sample from beta CDF
        index,beta = sampleCDF(CDF_beta_vals,betaCDF)
        # For this beta value, get our alpha CDF
        alphaCDF = getAlphaCDF(nAlpha,len(fullBetas),index,fullSAB,alphas)
        
        # If the valid alpha values for this chosen beta are not going to occur
        # we should just continue
        aMin, aMax = getAlphaMinMaxIndices(E,beta,kb,T,A,alphas)
        if (alphaCDF[aMax] < 1e-12): continue

        # Subtrack out the non-physical values
        H = [None]*nAlpha
        for a in range(nAlpha):
            if   (alphas[a] <= alphas[aMin]): H[a] = 0
            elif (alphas[a] >= alphas[aMax]): H[a] = 1
            else:  H[a] = ( (alphaCDF[a]-alphaCDF[aMin]) / \
                            (alphaCDF[aMax]-alphaCDF[aMin]) )
        index,alpha = sampleCDF(alphas,H)

        E_out_vec.append((beta*kb*T+E))

    plt.hist(E_out_vec,density=True,bins=np.logspace(-4.5,0,80),alpha=0.3)

"""
plt.gca().set_xscale("log")
plt.yscale('log', nonposy='clip')
axes = plt.gca()
plt.xlabel("E' [eV]")
plt.ylabel("Probability")

axes.set_xlim([4e-5,1e1])


plt.show()
"""












