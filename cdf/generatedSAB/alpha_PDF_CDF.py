import numpy as np
import matplotlib.pyplot as plt
from math import pi
from test09_sab import *
from freeGasSAB import *

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


def calcAlphaPDF(A,E,T,beta):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in test9_alphas if aMin <= x <= aMax]

    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = calcSym(alphas[a],beta)
        sabR = calcSym(alphas[a+1],beta)
        #sabL = getSABval(sab,a,63,n_beta)
        #sabR = getSABval(sab,a+1,63,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = calcSym(alphas[a],beta)/denominator
        #eq15[a] = getSABval(sab,a,63,n_beta)/denominator
    return alphas,eq15

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16




def PDF_CDF_at_various_temps(A,E,temps,beta):
    alpha_vecs, eq15_vecs, eq16_vecs = [], [], []
    for T in temps:
        alphas, eq15 = calcAlphaPDF(A=A,E=E,T=T,beta=beta)
        eq16 = calcAlphaCDF(eq15,alphas)
        alpha_vecs.append(alphas)
        eq15_vecs.append(eq15)
        eq16_vecs.append(eq16)

    for i in range(len(temps)):
        plt.plot(alpha_vecs[i],eq15_vecs[i])
    plt.show()

    for i in range(len(temps)):
        plt.plot(alpha_vecs[i],eq16_vecs[i])
    plt.show()

A = 18.0
#cA = 1.0
E = 8.0 
beta = 10
#temps = [300.0,650.0,1000.0,1500.0,2000.0] # Graphite
temps = [300.0,475.0,650.0,825.0,1000.0]   # Water
PDF_CDF_at_various_temps(A=A,E=E,temps=temps,beta=beta)






