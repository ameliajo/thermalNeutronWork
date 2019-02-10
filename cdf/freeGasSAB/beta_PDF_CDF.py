import numpy as np
import matplotlib.pyplot as plt
from math import pi


alphas = [.01008, .015, .0252, .033, 0.050406, .0756, 0.100812, 0.151218, 
  0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 0.504060, 
  0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
  1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 1.873790, 
  2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 3.578320, 
  3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 6.816500, 
  7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 12.96740, 
  14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 24.65260, 
  27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 46.85030, 50.0]

betas = [0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, 
  .0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 0.564547, 
  0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 1.129090, 1.209740, 
  1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 1.693640, 1.774290, 1.854940, 
  1.935590, 2.016240, 2.096890, 2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 
  2.580790, 2.669500, 2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 
  3.595360, 3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 
  5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 
  10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 17.17330, 18.72180, 
  20.42450, 22.29760, 24.35720, 25.0]




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
            sabL = getSABval(freeSAB,a,b,len(fullBetas))
            sabR = getSABval(freeSAB,a+1,b,len(fullBetas))
            #sabL = calcAsym(alphas[a],fullBetas[b])
            #sabR = calcAsym(alphas[a+1],fullBetas[b])
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



n_beta = len(betas)
freeSAB = [0.0]*(2*len(betas)-1)*len(alphas)
fullBetas = [-x for x in betas[1:]][::-1] + betas
for a in range(len(alphas)):
    for b in range(len(fullBetas)):
        freeSAB[a*len(fullBetas)+b] = (calcAsym(alphas[a],fullBetas[b]))




A = 18.0
E = 1.0 
temps = [300.0,475.0,650.0,825.0,1000.0]   # Water
PDF_CDF_at_various_temperatures(A,E,temps,fullBetas)
