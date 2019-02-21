import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *

continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]

import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('beta values')
    #plt.title(title)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    #plt.yscale('log')
    plt.show()



#alphas = [.01008, .015, .0252, .033, 0.050406, .0756, 0.100812, 0.151218, 
#  0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 0.504060, 
#  0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
#  1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 1.873790, 
#  2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 3.578320, 
#  3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 6.816500, 
#  7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 12.96740, 
#  14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 24.65260, 
#  27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 46.85030, 50.0]

#betas = [0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, 
#  .0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 0.564547, 
#  0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 1.129090, 1.209740, 
#  1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 1.693640, 1.774290, 1.854940, 
#  1.935590, 2.016240, 2.096890, 2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 
#  2.580790, 2.669500, 2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 
#  3.595360, 3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 
#  5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 
#  10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 17.17330, 18.72180, 
#  20.42450, 22.29760, 24.35720, 25.0]

#alphas = list(np.linspace(0.0001,60,100))
#betas = list(np.linspace(0.0,40,100))
alphas = list(np.linspace(0.01,60,500))
betas = list(np.linspace(0.0,40,100))
#alphas = list(np.linspace(1,30,3))
#betas = list(np.linspace(0.0,40,4))




n_alpha = len(alphas)
n_beta = len(betas)

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]

fullRedo = True
fullRedo = False
width = None 
temps = [296.0]   
temps = [296.0,475.0,650.0,825.0,1000.0]   

temps = list(np.linspace(296,1000,20))

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



def getFullSAB(alphas,betas,temps,fullBetas,sabs): 
    #freeSAB = [[0.0]*(2*len(betas)-1)*len(alphas) for t in range(len(temps))]
    fullSAB = [[0.0]*(2*len(betas)-1)*len(alphas) for t in range(len(temps))]
    for t in range(len(temps)):
        for a in range(len(alphas)):
            for b in range(len(fullBetas)):
                b_for_positive_betas = abs(b-len(betas)+1)
                assert(abs(fullBetas[b]) == betas[b_for_positive_betas])
                # This creates the symmetric SAB
                #freeSAB[t][a*len(fullBetas)+b] = calcSym(alphas[a],fullBetas[b])
                fullSAB[t][a*len(fullBetas)+b] = getSABval(sabs[t],a,b_for_positive_betas,len(betas))

                # We make it asymmetric to make sure it matches eq14
                #freeSAB[t][a*len(fullBetas)+b] *= np.exp(-fullBetas[b]*0.5)
                if (fullBetas[b] > 0): fullSAB[t][a*len(fullBetas)+b] *= np.exp(-fullBetas[b])
    return fullSAB



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



def calcIntegralAcrossAlpha(A,E,t,fullSAB,fullBetas,b):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,fullBetas[b],kb,temps[t],A)
    denominator = 0.0
    for a in range(len(alphas)-1):
        if aMin <= alphas[a] <= aMax:
            sabL = getSABval(fullSAB[t],a,b,len(fullBetas))
            sabR = getSABval(fullSAB[t],a+1,b,len(fullBetas))
            denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    return denominator

def calcBetaPDF(A,E,t,fullSAB,fullBetas,bMin,bMax):

    denominator = 0.0
    delta0 = fullBetas[bMin+1]-fullBetas[bMin]
    deltaN = fullBetas[bMax-1]-fullBetas[bMax-2]
    denominator = calcIntegralAcrossAlpha(A,E,t,fullSAB,fullBetas,bMin)*delta0*0.5 +\
                  calcIntegralAcrossAlpha(A,E,t,fullSAB,fullBetas,bMax-1)*deltaN*0.5
    for b in range(bMin+1,bMax-1):
        deltaL = fullBetas[b]-fullBetas[b-1]
        deltaR = fullBetas[b+1]-fullBetas[b]
        denominator += calcIntegralAcrossAlpha(A,E,t,fullSAB,fullBetas,b) * (deltaL+deltaR) * 0.5

    eq14 = [ calcIntegralAcrossAlpha(A,E,t,fullSAB,fullBetas,b) / denominator \
           for b in range(bMin,bMax-1)]
    return eq14



def calcBetaCDF(betas,eq14,bMin):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5*(fullBetas[bMin+b]-fullBetas[bMin+b-1]) 
    return eq16


def PDF_CDF_at_various_temperatures(A,E,temps,fullSAB,fullBetas):
    betas_vec = []
    eq14_vec = []
    eq16_vec = []
    for t in range(len(temps)):
        bMin,bMax = getValidBetasRange(A,E,temps[t],fullBetas)
        eq14 = calcBetaPDF(A,E,t,fullSAB,fullBetas,bMin,bMax)
        eq16 = calcBetaCDF(fullBetas,eq14,bMin)
        betas_vec.append(fullBetas[bMin:bMax-1])
        eq14_vec.append(eq14)
        eq16_vec.append(eq16)
    return betas_vec,eq16_vec


A = 1.0
E = 1.0 

fullBetas = [-x for x in betas[1:]][::-1] + betas
sabs_NJOY = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=True,fullRedo=fullRedo,width=width,oscE=oscE,oscW=oscW) for T in temps]
fullSAB = getFullSAB(alphas,betas,temps,fullBetas,sabs_NJOY)
beta_vecs,cdf = PDF_CDF_at_various_temperatures(A,E,temps,fullSAB,fullBetas)

betasToLookFor = list(np.linspace(-40,18,20))
scalarMap, colorBar = prepPlot(betasToLookFor)
for b,betaToLookFor in enumerate(betasToLookFor):
    tempVals = []
    cdfVals = []
    for t in range(len(temps)):
        bPrime = None 
        for i in range(len(beta_vecs[t])-1):
            if beta_vecs[t][i] <= betaToLookFor <= beta_vecs[t][i+1]: bPrime = i
        if bPrime != None:
            tempVals.append(temps[t])
            cdfVals.append(cdf[t][bPrime])
    plt.plot(tempVals,cdfVals,color=scalarMap.to_rgba(b))
    plt.plot(tempVals,cdfVals,'x',markersize=2,color=scalarMap.to_rgba(b))



plt.xlabel("Temperature (K)")
plt.ylabel("CDF")
finishPlotting(colorBar,"")


#sabs_MINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=False,fullRedo=False,width=width,oscE=oscE,oscW=oscW) for T in temps]
#fullSAB = getFullSAB(alphas,betas,temps,fullBetas,sabs_MINE)
#PDF_CDF_at_various_temperatures(A,E,temps,fullSAB,fullBetas)
#plt.show()


