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

alphas = list(np.linspace(0.01,60,500))
betas = list(np.linspace(0.0,40,100))

n_alpha = len(alphas)
n_beta = len(betas)

kb = 8.6173303e-5
oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]

def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax

def getAlphaRange(E,beta,T,A,alphas):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in alphas if aMin <= x <= aMax]
    return alphas

def getAlphaMinMaxIndices(E,beta,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    aMin_index = 0
    aMax_index = 0
    for a in range(len(alphas)):
        if alphas[a] > aMin and aMin_index == 0:
            aMin_index = a
        if alphas[a] <= aMax:
            aMax_index = a
    return aMin_index,aMax_index

def calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,index,n_beta)
        sabR = getSABval(sabs[t],a+1,index,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = getSABval(sabs[t],a,index,n_beta)/denominator
    return eq15

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16

def PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabs):
    eq19_vecs = []
    for t in range(len(temps)):
        eq18 = calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs)
        eq19 = calcAlphaCDF(eq18,alphas)
        eq19_vecs.append(eq19)
    return eq19_vecs

fullRedo = False
fullRedo = True
runNJOY = False
width = None 
temps = [296.0,475.0,650.0,825.0,1000.0]   # Water



A = 0.98
E = 1.0 

index = 25
beta = betas[index]
print("BETA",beta)
fullRedo = True
fullRedo = False

sabsMINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=True,fullRedo=fullRedo,\
            width=width,oscE=oscE,oscW=oscW) for T in temps]
eq19_vecs = PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabsMINE)



for t in range(len(temps)):
    aMin_index, aMax_index = getAlphaMinMaxIndices(E,beta,temps[t],A,alphas)
    alphasToPlot = []
    H = []
    H_hat = []
    for a in range(aMin_index,aMax_index+1):
        alphasToPlot.append(alphas[a])
        H.append( (eq19_vecs[t][a]-eq19_vecs[t][aMin_index]) / \
                  (eq19_vecs[t][aMax_index]-eq19_vecs[t][aMin_index]) )

    plt.plot(alphasToPlot,H)

plt.xlabel("alpha")
plt.ylabel("CDF")
plt.show()









































