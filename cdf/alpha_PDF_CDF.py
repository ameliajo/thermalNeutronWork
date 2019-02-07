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


def plotAlphaPDF(A,E,T,beta):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = list(np.linspace(aMin,aMax,1000))
    print(T,aMin,aMax)

    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = calcSym(alphas[a],beta)
        sabR = calcSym(alphas[a+1],beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = calcSym(alphas[a],beta)/denominator

    plt.plot(alphas,eq15)
    #print(np.trapz(eq15, x=alphas))


A = 18.0 
E = 7.0    # 1 eV
beta = 10
plotAlphaPDF(A=A,E=E,T=300.0,beta=beta)
plotAlphaPDF(A=A,E=E,T=475.0,beta=beta)
plotAlphaPDF(A=A,E=E,T=650.0,beta=beta)
plotAlphaPDF(A=A,E=E,T=825.0,beta=beta)
plotAlphaPDF(A=A,E=E,T=1000.0,beta=beta)
plt.show()


