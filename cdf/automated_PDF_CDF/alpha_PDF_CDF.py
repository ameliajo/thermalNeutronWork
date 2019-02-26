import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *
from rho import continRho

n_alpha, n_beta = 500, 100
alphas = list(np.linspace(0.01,60,n_alpha))
betas = list(np.linspace(0.0,40,n_beta))

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]


def getSABval(sab,a,b,n_beta): return sab[a*n_beta+b]

def getAlphaMinMax(E,beta,T,A):
    kb = 8.6173303e-5
    return ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T ), \
           ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )

def getAlphaRange(E,beta,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,T,A)
    return [x for x in alphas if aMin <= x <= aMax]

def calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,index,n_beta)
        sabR = getSABval(sabs[t],a+1,index,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    return [getSABval(sabs[t],a,index,n_beta)/denominator for a in range(len(alphas))]

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16

def PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabs):
    alpha_vecs, eq15_vecs, eq16_vecs = [], [], []
    for t in range(len(temps)):
        theseAlphas = getAlphaRange(E,beta,temps[t],A,alphas)
        eq15 = calcAlphaPDF(A,E,t,beta,theseAlphas,n_beta,index,sabs)
        alpha_vecs.append(theseAlphas)
        eq15_vecs.append(eq15)
        eq16_vecs.append(calcAlphaCDF(eq15,theseAlphas))

    for i in range(len(temps)): 
        plt.plot(alpha_vecs[i],eq15_vecs[i],label=str(temps[i])+' K')
    plt.legend(loc='best')
    plt.xlabel('alpha')
    plt.ylabel('PDF')
    plt.show()

    for i in range(len(temps)): 
        plt.plot(alpha_vecs[i],eq16_vecs[i],label=str(temps[i])+' K')
    plt.legend(loc='best')
    plt.xlabel('alpha')
    plt.ylabel('CDF')
    plt.show()



A = 0.98
E = 1.0 

# A with 18 and E with 10 looks great
runNJOY = False
width = None 
temps = [296.0,475.0,650.0,825.0,1000.0]   # Water


index = 25
beta = betas[index]
print("BETA",beta)
fullRedo = False

sabsMINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=runNJOY,fullRedo=False,\
            width=width,oscE=oscE,oscW=oscW) for T in temps]
PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabsMINE)





















