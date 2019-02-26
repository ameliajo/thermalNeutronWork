import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *
from rho import continRho


n_alpha, n_beta = 500, 100

alphas = list(np.linspace(0.01,60,n_alpha))
betas  = list(np.linspace(0.0,40,n_beta))

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]

def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def getAlphaMinMax(E,beta,T,A):
    kb = 8.6173303e-5
    return ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T ), \
           ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )


def getAlphaMinMaxIndices(E,beta,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,T,A)
    aMin_index, aMax_index = 0, 0
    for a in range(len(alphas)):
        if alphas[a] > aMin and aMin_index == 0: aMin_index = a
        if alphas[a] <= aMax:                    aMax_index = a
    return aMin_index,aMax_index

def calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,index,n_beta)
        sabR = getSABval(sabs[t],a+1,index,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])
    return [getSABval(sabs[t],a,index,n_beta)/denominator for a in range(len(alphas))]

def calcAlphaCDF(eq18,alphas):
    eq19 = [0.0]*len(eq18)
    for a in range(1,len(alphas)):
        eq19[a] = eq19[a-1]+(eq18[a]+eq18[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq19
    

def PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabs):
    eq19_vecs = []
    for t in range(len(temps)):
        eq18 = calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs)
        eq19_vecs.append(calcAlphaCDF(eq18,alphas))
    return eq19_vecs


width = None 
temps = [296.0,475.0,650.0,825.0,1000.0]   # Water
A = 0.98
E = 1.0 
index = 25
beta = betas[index]
print("BETA",beta)

sabsMINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=True,fullRedo=False,\
            width=width,oscE=oscE,oscW=oscW) for T in temps]
eq19_vecs = PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabsMINE)



for t in range(len(temps)):
    aMin_index, aMax_index = getAlphaMinMaxIndices(E,beta,temps[t],A,alphas)
    alphasToPlot, H = [], []
    for a in range(aMin_index,aMax_index+1):
        alphasToPlot.append(alphas[a])
        H.append( (eq19_vecs[t][a]-eq19_vecs[t][aMin_index]) / \
                  (eq19_vecs[t][aMax_index]-eq19_vecs[t][aMin_index]) )
    plt.plot(alphasToPlot,H,label=str(temps[t])+' K')

plt.legend(loc='best')
plt.xlabel("alpha")
plt.ylabel("CDF")
plt.show()









































