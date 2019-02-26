from numpy import exp
import numpy as np
import matplotlib.pyplot as plt



def interpolate(xList,yList,x):
    for i in range(len(xList)-1):
        if (xList[i] <= x < xList[i+1]):
            m = (yList[i+1]-yList[i]) / (xList[i+1]-xList[i])
            b = yList[i] - m*xList[i]
            return m*x+b
    return 0.0 if (xList[-1] < x) or xList[0] > x else yList[-1]

def getVal(betaGrid,Tn,b):
    if abs(b) >= len(betaGrid):
        return 0.0
    if b < 0:
        return Tn[abs(b)]*exp(betaGrid[abs(b)])
    return Tn[b]



inputBetas = list(np.linspace(0.0,4,10)) + list(np.linspace(1.0,3,25))
inputBetas.sort()
inputT1 = [-2*x*x+8*x for x in inputBetas]

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


# Lets assume an alpha = 0.05 and a lambda = 1.2
beta = np.linspace(0,10,21) 
N = len(beta)
T1 = [interpolate(inputBetas,inputT1,abs(beta[i])) for i in range(N)]
Tn = T1[:]
S = [0.0]*N
S[0] = exp(-0.05*1.2)
for n in range(1,15):
    alpha_n_term = (exp(-0.05*1.2)*(0.05*1.2)**n)*(1.0/np.math.factorial(n))
    for b in range(N):
        S[b] = alpha_n_term * Tn[b]
    Tn = convolve(beta,T1,Tn)
    plt.plot(beta,S)


plt.show()
























