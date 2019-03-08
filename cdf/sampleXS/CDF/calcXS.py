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
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
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
            val = frac*(x[i+1]-x[i])+x[i]
            return i,val



def getAlphaCDF(nAlpha,nFullBetas,index,fullSAB,alphas,aMin):
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

nAlpha, nBeta = 100, 1000
alphas = list(np.linspace(0.0001,30,nAlpha))
betas  = list(np.linspace(0.00,300,nBeta))
betas  = list(np.logspace(-6,2.1,nBeta))

useNJOY  = True
fullRedo = False  
fullRedo = True
width    = None 


fullBetas = [-x for x in betas[1:]][::-1] + betas
SAB = getSAB(alphas,betas,T,continRho,useNJOY,fullRedo,width,oscE,oscW) 
fullSAB = getFullSAB(alphas,betas,fullBetas,SAB)

N = 500
N = int(100000)
N = 50000

Energies = [0.0253]
Energies = [0.005,0.01,0.95]
Energies = [0.01]
Energies = [3.0]
Energies = [0.005]
Energies = [0.0253]
Energies = [3.12]
Energies = [0.0005,0.0253,0.2907,0.95,3.12]
#scalarMap, colorBar = prepPlot(Energies)


for i,E in enumerate(Energies):
    E_out_vec = []
    betas_vec= []
    print("E = ",E)

    relevantBetas, betaCDF = PDF_CDF_at_various_temperatures(A,E,T,fullSAB,fullBetas)

    xsi_Vals = []
    for n in range(N):
        if n in [1000,3000,10000,50000]:
            print(n)
        index,beta = sampleCDF(relevantBetas,betaCDF)

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
        #betas_vec.append((beta))


    #sns.kdeplot(E_out_vec,shade=True);
    #print(min(betas_vec),max(betas_vec))
    #print(min(E_out_vec),max(E_out_vec))
    #plt.hist(E_out_vec,bins=np.logspace(-3,1,500),normed=True,alpha=0.3)
    #sns.kdeplot(E_out_vec,label=str(E)+' eV',shade=True,cut=10);
    #sns.distplot(E_out_vec,bins=np.logspace(-3,1,2),hist=True,label=str(E)+' eV');
    plt.hist(E_out_vec,normed=True,bins=np.logspace(-3,1,100),alpha=0.3)
    #plt.hist(E_out_vec,bins=np.logspace(-3,1,500),normed=True,alpha=0.3)
    #plt.hist(E_out_vec,bins=np.linspace(-40,5,100),alpha=0.3)
    #plt.hist(E_out_vec,bins=np.linspace(0.0,3.5,100),alpha=0.3)

#sns.kdeplot(E_out_vec);
plt.gca().set_xscale("log")
plt.yscale('log', nonposy='clip')

#plt.xscale('log')
#plt.yscale('log')

axes = plt.gca()
axes.set_xlim([1e-3,1e1])
axes.set_ylim([1e-3,1e3])


plt.show()




#plt.plot(relevantBetas,betaCDF)
#plt.xlabel('alpha')
#plt.ylabel('beta')

#finishPlotting(colorBar,"")











