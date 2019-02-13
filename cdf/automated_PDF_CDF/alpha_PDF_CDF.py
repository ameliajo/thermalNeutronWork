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

alphas = [.01008, .015, .0252, .033, 0.050406, .0756, 0.100812, 0.151218, 
  0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 0.504060, 
  0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 
  1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 1.873790, 
  2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 3.578320, 
  3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 6.816500, 
  7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 12.96740, 
  14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 24.65260, 
  27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 46.85030, 50.0]

betas = [0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, 
  .0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 0.564547, 
  0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 1.129090, 1.209740, 
  1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 1.693640, 1.774290, 1.854940, 
  1.935590, 2.016240, 2.096890, 2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 
  2.580790, 2.669500, 2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 
  3.595360, 3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 
  5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 
  10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 17.17330, 18.72180, 
  20.42450, 22.29760, 24.35720, 25.0]

n_alpha = len(alphas)
n_beta = len(betas)

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]

fullRedo = False
fullRedo = True
width = None 
temps = [296.0,475.0,650.0,825.0,1000.0]   # Water

sabs = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=False,fullRedo=fullRedo,width=width,oscE=oscE,oscW=oscW) for T in temps]

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

def getAlphaRange(E,beta,T,A,alphas):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in alphas if aMin <= x <= aMax]
    return alphas

def calcAlphaPDF(A,E,t,beta,alphas,n_beta):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,21,n_beta)
        sabR = getSABval(sabs[t],a+1,21,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = getSABval(sabs[t],a,21,n_beta)/denominator
    return eq15

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16

def PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta):
    alpha_vecs, eq15_vecs, eq16_vecs = [], [], []
    for t in range(len(temps)):
        alphas = getAlphaRange(E,beta,temps[t],A,alphas)
        eq15 = calcAlphaPDF(A=A,E=E,t=t,beta=beta,alphas=alphas,n_beta=n_beta)
        alpha_vecs.append(alphas)
        eq15_vecs.append(eq15)
        eq16 = calcAlphaCDF(eq15,alphas)
        eq16_vecs.append(eq16)

    for i in range(len(temps)): plt.plot(alpha_vecs[i],eq15_vecs[i])
    plt.show()

    for i in range(len(temps)): plt.plot(alpha_vecs[i],eq16_vecs[i])
    plt.show()



A = 1.0
E = 1.0 

# A with 18 and E with 10 looks great
beta = 10
beta = betas[21]
PDF_CDF_at_various_temps(A=A,E=E,temps=temps,beta=beta,alphas=alphas,n_beta=n_beta)



