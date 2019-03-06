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
        Tn = convolve(beta,T1,Tn,len(beta))
        n += 1
        if (len(diffs)>5 and diffs[-1] < 1e-6):
            return diffs
    return diffs

    #for n in range(1,numIter):
    #    aTerms *= alpha / (1.0*n)
    #    diffs[n-1] = sum([aTerms*Tn[b]*(beta[b+1]-beta[b]) \
    #                      for b in range(len(beta)-1)])
    #    Tn = convolve(beta,T1,Tn,len(beta))
    #return diffs

numIter = 20
numBeta = 300
beta = np.linspace(0,20,numBeta) 


alphas = [0.1,0.2,0.5,1.0,2.0,2.5,3.0]
alphas = [0.1,0.5,1.0,2.5,4.0,5.0,10.0]

peakLocations = []
convergenceLocations1 = []
convergenceLocations2 = []
convergenceLocations3 = []
convergenceLocations4 = []
convergenceLocations5 = []

def getConvergenceLocation(location,diffs,tol):
    for i in range(location,len(diffs)):
        if diffs[i] < tol:
            return i

tols = [1e-6,1e-4,1e-3,1e-2,1e-1]
tols = [1e-6,20,40,80,120]
for a,alpha in enumerate(alphas):
    diffs = getConvergenceOfSum(alpha,beta,numIter,rhoBeta)
    #plt.plot([x/sum(diffs) for x in diffs],label=str(alpha))
    peakLocations.append(diffs.index(max(diffs)))
    peak = peakLocations[-1]
    convergenceLocations1.append(len(diffs))
    convergenceLocations2.append(getConvergenceLocation(peak,diffs,tols[1]))
    convergenceLocations3.append(getConvergenceLocation(peak,diffs,tols[2]))
    convergenceLocations4.append(getConvergenceLocation(peak,diffs,tols[3]))
    convergenceLocations5.append(getConvergenceLocation(peak,diffs,tols[4]))







colors = ['#ff0000', '#ff4000', '#ff8000', '#ffbf00', '#ffff00']


#plt.xlabel('# Iterations')
#plt.ylabel('Contribution to S(a,b)')
#plt.legend(loc='best')
#plt.show()
#plt.plot(alphas,convergenceLocations1E6,'bo')
plt.plot(alphas,convergenceLocations1,colors[0],label=str(tols[0])+'%')
plt.plot(alphas,convergenceLocations2,colors[1],label=str(tols[1])+'%')
plt.plot(alphas,convergenceLocations3,colors[2],label=str(tols[2])+'%')
plt.plot(alphas,convergenceLocations4,colors[3],label=str(tols[3])+'%')
plt.plot(alphas,convergenceLocations5,colors[4],label=str(tols[4])+'%')
#plt.plot(alphas,peakLocations,'k',label='peak')
#plt.plot(alphas,convergenceLocations1E3,'go')
#plt.plot(alphas,convergenceLocations1E3,'g',label='convergence (1e-3)')
plt.legend(loc='best')
plt.xlabel('alpha')
#plt.ylabel('# Iterations necessary to get to Peak')
plt.ylabel('# Iterations')
plt.show()
#plt.plot(convergenceLocations)
#plt.show()
"""
"""























"""


numIter = 50
scalarMap, colorBar = prepPlot(list(range(numIter)))

alpha = 3.0
numBeta = 100
beta = np.linspace(0,3,numBeta) 

diffs = getConvergenceOfSum(alpha,beta,numIter)





trueDiffs = [2.4121836009801084e-07, 2.2063872462386213e-06, 1.4594017441093847e-05, 7.4082832786500933e-05, 0.00030335272399319212, 0.0010392399199991198, 0.0030588047305379882, 0.0078899959373381675, 0.01811098956690866, 0.037447640715426382, 0.070438204669199797, 0.12151807124190461, 0.19360065951640482, 0.28651800519207349, 0.39588864651986738, 0.51296215871433704, 0.62570981687186966, 0.72098897017375563, 0.78719854279971324, 0.81664985299723114, 0.80698068786804744, 0.76128216839348695, 0.68703088370327159, 0.59425387602170132, 0.49349585949127134, 0.39409651135804924, 0.30308820384770885, 0.22478959687472144, 0.16098148749315844, 0.11145067692707357, 0.074675366176478616, 0.048474155171743728, 0.030514346663416101, 0.018644741008298982, 0.011067297813397185, 0.0063872401981762581, 0.0035867856509832724, 0.0019612567202264512, 0.0010449613497483085, 0.00054285960309393539, 0.0002751484361551681, 0.0001361433657559043, 6.579931947040709e-05, 3.1079648883037009e-05, 1.4354377732747896e-05, 6.485749329879233e-06, 2.8681921131320359e-06, 1.242009045359222e-06, 5.268631199836487e-07]
for i in range(len(diffs)):
    assert(abs((diffs[i]-trueDiffs[i])/trueDiffs[i])<1e-14)


print(max(diffs))
"""












