import numpy as np
import matplotlib.pyplot as plt
from math import pi



test9_alphas = [.01008, .015, .0252, .033, 0.050406, .0756, 0.100812, 0.151218, 
  0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 0.504060, 
  0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
  1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 1.873790, 
  2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 3.578320, 
  3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 6.816500, 
  7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 12.96740, 
  14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 24.65260, 
  27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 46.85030, 50.0]


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

def getAlphaRange(E,beta,T,A):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in test9_alphas if aMin <= x <= aMax]
    return alphas

def calcAlphaPDF(A,E,T,beta,alphas):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = calcSym(alphas[a],beta)
        sabR = calcSym(alphas[a+1],beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = calcSym(alphas[a],beta)/denominator
    return eq15

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16

def PDF_CDF_at_various_temps(A,E,temps,beta):
    alpha_vecs, eq15_vecs, eq16_vecs = [], [], []
    for T in temps:
        alphas = getAlphaRange(E,beta,T,A)
        eq15 = calcAlphaPDF(A=A,E=E,T=T,beta=beta,alphas=alphas)
        alpha_vecs.append(alphas)
        eq15_vecs.append(eq15)
        eq16 = calcAlphaCDF(eq15,alphas)
        eq16_vecs.append(eq16)

    for i in range(len(temps)): plt.plot(alpha_vecs[i],eq15_vecs[i])
    plt.show()

    for i in range(len(temps)): plt.plot(alpha_vecs[i],eq16_vecs[i])
    plt.show()

A = 18.0
E = 10.0 
beta = 10
#temps = [300.0,650.0,1000.0,1500.0,2000.0] # Graphite
temps = [300.0,475.0,650.0,825.0,1000.0]   # Water
PDF_CDF_at_various_temps(A=A,E=E,temps=temps,beta=beta)






