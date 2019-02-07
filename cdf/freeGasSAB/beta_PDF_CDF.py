import numpy as np
import matplotlib.pyplot as plt
from math import pi


def calcAsym(alpha,beta):
    return (1.0/(4.0*pi*alpha)**0.5) * np.exp( - (alpha+beta)**2 / (4.0*alpha) )


def calcSym(alpha,beta):
    return np.exp(beta*0.5)*calcAsym(alpha,beta)

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax


def calcIntegralAcrossAlpha(A,E,T,beta):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = list(np.linspace(aMin,aMax,50))
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = calcAsym(alphas[a],beta)
        sabR = calcAsym(alphas[a+1],beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    return denominator

def calcBetaPDF(A,E,T):
    kb = 8.6173303e-5
    bMin = -E/(kb*T)
    bMax = 20.0
    betas = list(np.linspace(bMin,bMax,100))

    denominator = 0.0
    for b in range(len(betas)-1):
        denominator += ( calcIntegralAcrossAlpha(A,E,T,betas[b]) +    \
                         calcIntegralAcrossAlpha(A,E,T,betas[b+1]) ) *\
                       0.5 * ( betas[b+1]-betas[b] )
    eq14 = [0.0]*len(betas)
    for b in range(len(betas)):
        eq14[b] = calcIntegralAcrossAlpha(A,E,T,betas[b]) / denominator
    return betas, eq14

def calcBetaCDF(betas,eq14):
    eq16 = [0.0]*len(betas)
    for b in range(1,len(betas)):
        eq16[b] = eq16[b-1] + (eq14[b-1]+eq14[b])*0.5*(betas[b]-betas[b-1]) 
    return eq16


def PDF_CDF_at_various_temperatures(A,E,temps):
    betas_vec = []
    eq14_vec = []
    eq16_vec = []
    for T in temps:
        betas, eq14 = calcBetaPDF(A,E,T)
        eq16 = calcBetaCDF(betas,eq14)
        betas_vec.append(betas)
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

