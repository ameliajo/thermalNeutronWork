import numpy as np
import matplotlib.pyplot as plt
from math import pi
from test09_sab import *
#from freeGasSAB import *

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

freeSAB_pos = []
freeSAB_neg = []
n_alpha = len(alphas)
n_beta = len(betas)
freeSAB = [0.0]*(2*n_beta-1)*n_alpha #for x in range(n_alpha)]
fullBetas = [-x for x in betas[1:]][::-1] + betas
for a in range(n_alpha):
    for b in range(n_beta):
        freeSAB_pos.append(calcAsym(alphas[a],betas[b]))
        freeSAB_neg.append(calcAsym(alphas[a],-betas[b]))
for a in range(n_alpha):
    for b in range(2*n_beta-1):
        freeSAB[a*len(fullBetas)+b] = (calcAsym(alphas[a],fullBetas[b]))





def calcIntegralAcrossAlpha(A,E,T,beta,b):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    denominator = 0.0
    for a in range(len(alphas)-1):
        if aMin <= alphas[a] <= aMax:
            sabL = calcAsym(alphas[a],fullBetas[b])
            sabR = calcAsym(alphas[a+1],fullBetas[b])
            denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    return denominator

def calcBetaPDF(A,E,T):
    kb = 8.6173303e-5
    bMin = -E/(kb*T)
    bMax = 20.0
    denominator = 0.0
    betasToPlot = []
    for b in range(len(fullBetas)-1):
        if bMin <= fullBetas[b] <= bMax:
            denominator += ( calcIntegralAcrossAlpha(A,E,T,fullBetas[b],b) +    \
                             calcIntegralAcrossAlpha(A,E,T,fullBetas[b+1],b+1) ) *\
                           0.5 * ( fullBetas[b+1]-fullBetas[b] )

    eq14 = []
    for b in range(len(fullBetas)):
        if bMin <= fullBetas[b] <= bMax:
            betasToPlot.append(fullBetas[b])
            eq14.append(calcIntegralAcrossAlpha(A,E,T,fullBetas[b],b) / denominator)
    return betasToPlot,eq14

def calcBetaCDF(betas,eq14):
    eq16 = [0.0]*len(eq14)
    for b in range(1,len(eq14)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5*(fullBetas[b]-fullBetas[b-1]) 
    return eq16


def PDF_CDF_at_various_temperatures(A,E,temps):
    betas_vec = []
    eq14_vec = []
    eq16_vec = []
    for T in temps:
        betasToPlot, eq14 = calcBetaPDF(A,E,T)
        eq16 = calcBetaCDF(fullBetas,eq14)
        betas_vec.append(betasToPlot)
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
PDF_CDF_at_various_temperatures(A,E,temps)
