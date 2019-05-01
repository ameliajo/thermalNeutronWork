
import numpy as np
from math import factorial as f
import matplotlib.pyplot as plt
from scipy.special import iv



alpha = 1.0
betaMin = 0.0
betaMax = 10.0

beta_i = 0.5
betaGrid = np.linspace(betaMin,betaMax,(betaMax-betaMin)/beta_i,endpoint=False)

def P(beta,beta_i):
    if abs(beta-beta_i) < 1e-6:
        return 1.0/(2.0*beta_i*np.sinh(beta_i*0.5))
    return 0.0




lambda_s = np.cosh(beta_i*0.5) / (beta_i*np.sinh(beta_i*0.5))
time = np.linspace(-100,100,100)
fullBeta = list(betaGrid[1:][::-1])+ list(betaGrid)
print([P(beta,beta_i) for beta in fullBeta])

gamma = []
#for t in time:
#    gamma.append(alpha*np.trapz([P(beta,beta_i) * (1.0-np.exp(1j*beta*t)) \
#                                *np.exp(-beta/2) for beta in fullBeta],x=fullBeta))

#    print([P(beta,beta_i) for beta in fullBeta])


 

#plt.plot(gamma)
#plt.show()

print()
print()
print()



def getContinVersion(alpha,betaMin,betaMax,beta_i,betaGrid):
    time = np.linspace(-100,100,1000)
    lambda_s = np.cosh(beta_i*0.5) / (beta_i*np.sinh(beta_i*0.5))
    continSAB = [0.0]*len(betaGrid)
    for n,beta1 in enumerate(betaGrid):
        gamma = [alpha*np.trapz([P(beta,beta_i)*(1.0-np.exp(-1j*beta*t))*np.exp(-beta*0.5) for beta in betaGrid],x=betaGrid) for t in time]
        continSAB[n] = 1.0/(2*3.14159) * np.trapz([np.exp(1j*beta1*time[t])*np.exp(-gamma[t]) for t in range(len(time))],x=time)
    return continSAB



def getContinVersionGood(alpha,betaMin,betaMax,beta_i,betaGrid):
    lambda_s = np.cosh(beta_i*0.5) / (beta_i*np.sinh(beta_i*0.5))
    print(lambda_s)
    continSAB = [0.0]*len(betaGrid)
    for n,beta in enumerate(betaGrid):
        continSAB[n] = np.exp(-alpha*lambda_s-0.5*beta) * 1.0/f(n) * (alpha)**n * \
                       0.5 / (beta*np.sinh(beta_i*0.5))**n
    return continSAB





def getDiscreVersion(alpha,betaMin,betaMax,beta_i,betaGrid):
    besselArg = alpha / ( beta_i * np.sinh(beta_i*0.5) )
    lambda_s = 1.0/(np.tanh(beta_i*0.5) * beta_i)
    #print(lambda_s)
    discreSAB = [0.0]*len(betaGrid)
    for n,beta in enumerate(betaGrid):
        discreSAB[n] = np.exp(-alpha*lambda_s) * iv(n,besselArg) * np.exp(-n*beta_i*0.5)
    return discreSAB





#continSAB = getContinVersion(alpha,betaMin,betaMax,beta_i,betaGrid)
#continSABGood = getContinVersionGood(alpha,betaMin,betaMax,beta_i,betaGrid)
discreSAB = getDiscreVersion(alpha,betaMin,betaMax,beta_i,betaGrid)
lambda_s = np.cosh(beta_i*0.5) / (beta_i*np.sinh(beta_i*0.5))
print(np.exp(-alpha*lambda_s)*alpha*np.exp(beta_i*0.5)*P(beta_i,beta_i))
print(discreSAB[:3])
#plt.plot(betaGrid,continSAB,'o',label='contin')
#plt.plot(betaGrid,continSAB)
"""
plt.plot(betaGrid,continSABGood,'o',label='contin good')
plt.plot(betaGrid,continSABGood)
plt.plot(betaGrid,discreSAB,'o',label='discre')
plt.plot(betaGrid,discreSAB)
plt.legend(loc='best')
"""
"""
"""
plt.show()
#plt.plot(betaGrid,continSAB,'ro')
#plt.show()















