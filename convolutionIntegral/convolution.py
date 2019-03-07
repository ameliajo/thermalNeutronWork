from numpy import exp
import numpy as np
import matplotlib.pyplot as plt
from plot import *

continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]

rhoBeta = np.arange(0,0.00255*len(continRho),0.00255)


def getT1(continRho,alpha,rhoBeta):
    # Construct P(b) = rho(b) / 2*b*sinh(b/2)
    P  = [continRho[1]/rhoBeta[1]**2] 
    P += [continRho[b]/(2*rhoBeta[b]*np.sinh(rhoBeta[b]*0.5)) \
            for b in range(1,len(rhoBeta))]
    
    lambda_s = 0.0
    for i in range(len(rhoBeta)-1):
        lambda_s += ( P[i]*np.exp(-rhoBeta[i]*0.5) + \
                      P[i+1]*np.exp(-rhoBeta[i+1]*0.5) ) * \
                    ( rhoBeta[i+1]-rhoBeta[i] ) * 0.5
    
    inputT1 = [P[b]*exp(-rhoBeta[b]*0.5) for b in range(len(rhoBeta))] 
    return inputT1,lambda_s



def interpolate(xList,yList,x,lambda_s):
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


def convolve(beta,T1,Tlast,numBeta):
    Tn = [0.0]*numBeta
    for b in range(numBeta):
        for bp in range(numBeta-1):
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


def getConvergenceOfSum(alpha,beta,numIter,rhoBeta):
    inputT1,lambda_s = getT1(continRho,alpha,rhoBeta)
    beta = rhoBeta[:]
    T1 = [interpolate(rhoBeta,inputT1,abs(beta[i]),lambda_s) \
            for i in range(len(beta))]
    Tn = T1[:]
    diffs = [0.0]*(numIter-1)
    aTerms = exp(-alpha*lambda_s)

    sab = [0.0]*len(beta)
    n = 1
    diffs = []
    while True:
        aTerms *= alpha / (1.0*n)
        #diffs.append(sum([aTerms*Tn[b]*(beta[b+1]-beta[b]) \
        #                  for b in range(len(beta)-1)]))
        differences = 0.0
        for b in range(len(beta)):
            if (sab[b] > 1e-6 ):
                differences += aTerms*Tn[b] / sab[b]
            sab[b] = aTerms *Tn[b]
        diffs.append(differences)
        #diffs.append(aTerms*Tn[b])
        Tn = convolve(beta,T1,Tn,len(beta))
        n += 1
        if (len(diffs)>5 and diffs[-1] < 1e-4):
            return diffs
    return diffs

def getConvergenceLocation(location,diffs,tol):
    for i in range(location,len(diffs)):
        if diffs[i] < tol:
            return i



numIter = 20
numBeta = 300
beta = np.linspace(0,20,numBeta) 


alphas = [0.1,0.2,0.5,1.0,2.0,2.5,3.0]
alphas = [0.1,0.5,1.0,2.5,4.0,5.0,10.0]

plotConvergence = True
if plotConvergence:
    peakLocations = []
    convergenceLocations1 = []
    convergenceLocations2 = []
    convergenceLocations3 = []
    convergenceLocations4 = []
    convergenceLocations5 = []


    #tols = [1e-6,1e-4,1e-3,1e-2,1e-1]
    tols = [1e-6,20,40,80,120]
    convergences = [[None]*len(alphas) for x in range(len(tols))]
    for a,alpha in enumerate(alphas):
        diffs = getConvergenceOfSum(alpha,beta,numIter,rhoBeta)
        peakLocations.append(diffs.index(max(diffs)))
        peak = peakLocations[-1]
        convergenceLocations1.append(len(diffs))
        convergenceLocations2.append(getConvergenceLocation(peak,diffs,tols[1]))
        convergenceLocations3.append(getConvergenceLocation(peak,diffs,tols[2]))
        convergenceLocations4.append(getConvergenceLocation(peak,diffs,tols[3]))
        convergenceLocations5.append(getConvergenceLocation(peak,diffs,tols[4]))
        convergences[0][a] = len(diffs)
        convergences[1][a] = getConvergenceLocation(peak,diffs,tols[1])
        convergences[2][a] = getConvergenceLocation(peak,diffs,tols[2])
        convergences[3][a] = getConvergenceLocation(peak,diffs,tols[3])
        convergences[4][a] = getConvergenceLocation(peak,diffs,tols[4])



    colors = ['#ff0000', '#ff4000', '#ff8000', '#ffbf00', '#ffff00']


    plt.plot(alphas,convergences[0],colors[0],label=str(tols[0])+'%')
    plt.plot(alphas,convergences[1],colors[1],label=str(tols[1])+'%')
    plt.plot(alphas,convergences[2],colors[2],label=str(tols[2])+'%')
    plt.plot(alphas,convergences[3],colors[3],label=str(tols[3])+'%')
    plt.plot(alphas,convergences[4],colors[4],label=str(tols[4])+'%')
    plt.legend(loc='best')
    plt.xlabel('alpha')
    plt.ylabel('# Iterations')
    plt.show()




plotContributions = True
plotContributions = False
if plotContributions:
    for a,alpha in enumerate(alphas):
        diffs = getConvergenceOfSum(alpha,beta,numIter,rhoBeta)
        plt.plot([x/sum(diffs) for x in diffs],label='a = '+str(alpha))


    plt.xlabel('# Iterations')
    plt.ylabel('Contribution to S(a,b)')
    plt.legend(loc='best')
    plt.show()



