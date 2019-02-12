import sys
sys.path.append('NJOY_SAB_DATA')
import numpy as np
import matplotlib.pyplot as plt
from math import pi
#from test09_sab import *
from SAB_test09_296K import *
from SAB_test09_475K import *
from SAB_test09_650K import *
from SAB_test09_825K import *
from SAB_test09_1000K import *




def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def calcAsym(alpha,beta):
    return (1.0/(4.0*pi*alpha)**0.5) * np.exp( - (alpha+beta)**2 / (4.0*alpha) )

def calcSym(alpha,beta):
    return np.exp(beta*0.5)*calcAsym(alpha,beta)

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    print(aMin,aMax)
    return aMin,aMax
    #aMin = 0.003498
    #aMax = 88.22865
    #return 0.00349820, 88.228649

def getAlphaRange(E,beta,T,A,alphas):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in alphas if aMin <= x <= aMax]
    return alphas

def calcAlphaPDF(A,E,t,beta,alphas,n_beta):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,63,n_beta)
        sabR = getSABval(sabs[t],a+1,63,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = getSABval(sabs[t],a,63,n_beta)/denominator
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
    #plt.show()


temps = [296.0,475.0,650.0,825.0,1000.0]   # Water
assert(alphas_296 == alphas_475 == alphas_650 == alphas_825 == alphas_1000)
assert(betas_296 == betas_475 == betas_650 == betas_825 == betas_1000)

alphas = alphas_296[:]
betas = betas_296[:]
sabs = [sab_296,sab_475,sab_650,sab_825,sab_1000]
n_beta = len(betas)


A = 18.0
E = 10.0 

# A with 18 and E with 10 looks great
beta = 10
beta = betas[63]
PDF_CDF_at_various_temps(A=A,E=E,temps=temps,beta=beta,alphas=alphas,n_beta=n_beta)






