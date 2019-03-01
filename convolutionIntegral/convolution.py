from numpy import exp
import numpy as np
import matplotlib.pyplot as plt
from plot import *


# Lets assume an alpha value
alpha = 0.50
continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]

rhoBeta = np.arange(0,0.00255*len(continRho),0.00255)
P = [continRho[b]/(2*rhoBeta[b]*np.sinh(rhoBeta[b]*0.5)) for b in range(1,len(rhoBeta))]
P = [continRho[1]/rhoBeta[1]**2] + P

lambda_s = P[0] *2.0*np.cosh(rhoBeta[0 ]*0.5)*0.5*(rhoBeta[1] -rhoBeta[0] ) + \
           P[-1]*2.0*np.cosh(rhoBeta[-1]*0.5)*0.5*(rhoBeta[-1]-rhoBeta[-2])
for b in range(1,len(rhoBeta)-1):
    lambda_s += P[b]*2.0*np.cosh(rhoBeta[b]*0.5)*(rhoBeta[b]-rhoBeta[b-1])

inputT1 = [P[b]*exp(-rhoBeta[b]*0.5)/lambda_s for b in range(len(rhoBeta))] 

#plt.plot(rhoBeta,continRho)
#plt.plot(rhoBeta,P)
#plt.plot(rhoBeta,T1)
#plt.show()


def interpolate(xList,yList,x):
    for i in range(len(xList)-1):
        if (xList[i] <= x < xList[i+1]):
            m = (yList[i+1]-yList[i]) / (xList[i+1]-xList[i])
            b = yList[i] - m*xList[i]
            return m*x+b
    return 0.0 if (xList[-1] < x) or xList[0] > x else yList[-1]

def getVal(betaGrid,Tn,b):
    if abs(b) >= len(betaGrid): return 0.0
    if b < 0: return Tn[abs(b)]*exp(betaGrid[abs(b)])
    return Tn[b]


def convolve(beta,T1,Tlast):
    Tn = [0.0]*N
    for b in range(N):
        for bp in range(N-1):
            AL = T1[bp]*exp(beta[bp])
            AR = T1[bp+1]*exp(beta[bp+1])
            BL = getVal(beta,Tlast,b+bp)
            BR = getVal(beta,Tlast,b+bp+1)
            CL = T1[bp]
            CR = T1[bp+1]
            DL = getVal(beta,Tlast,b-bp)
            DR = getVal(beta,Tlast,b-bp-1)
            Tn[b] += ((AL*BL+CL*DL) + (AR*BR+CR*DR)) * \
                     (beta[bp+1]-beta[bp]) * 0.5
    return Tn


numIter = 10
scalarMap, colorBar = prepPlot(list(range(numIter)))
beta = np.linspace(0,3,300) 
N = len(beta)
T1 = [interpolate(rhoBeta,inputT1,abs(beta[i])) for i in range(N)]
Tn = T1[:]
S = [0.0]*N
S[0] = exp(-alpha*lambda_s)
diffs = []
for n in range(1,numIter):
    alpha_n_term = (exp(-alpha*lambda_s)*(alpha*lambda_s)**n) * \
                   (1.0/np.math.factorial(n))
    totalDiff = 0.0
    for b in range(N):
        S[b] += alpha_n_term * Tn[b]
        totalDiff += alpha_n_term * Tn[b] * 0.0001
    diffs.append(totalDiff)
    Tn = convolve(beta,T1,Tn)
    plt.plot(beta,S,color=scalarMap.to_rgba(n))
    if totalDiff < 1e-6:
        break
finishPlotting(colorBar)
#plt.plot(diffs)
#plt.show()
















