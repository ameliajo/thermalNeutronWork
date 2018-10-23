import matplotlib.pyplot as plt
from results_test_09_300K import *
#from results_test_09_300K_with_negative_beta import *
#from results_symBeta import *
from math import exp



def readSab(a,b,sab_water,betas):
    print("beta Value:",betas[b])
    print("sab Value: ",sab_water[a*len(betas)+b])
    print("sab Value: ",exp(betas[b]/2)*sab_water[a*len(betas)+b])





def getEq14(beta,E,T,Sab,A,alphaVec,lenBeta,b):
    alpha_min = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    alpha_max = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    #print(beta,alpha_min,alpha_max)
    validAlphaIndexes = [ i for i in range(len(alphaVec)) if \
                   (alphaVec[i] >= alpha_min and alphaVec[i] <= alpha_max) ]

    numerator = exp(-beta/2) * sum([Sab[i*lenBeta+b] for i in validAlphaIndexes])
    #numerator = sum([Sab[i*lenBeta+b] for i in validAlphaIndexes])
    #if beta < -24.5:
    #    print(exp(-beta/2),numerator)


    return numerator 



A = 18.02   # H2O = 1.01 + 1.01 + 16.00 = 18.02 
A = 1.0   # H2O = 1.01 + 1.01 + 16.00 = 18.02 
E = 1       # 1 eVc
T = 300     # 300 K
kb = 8.6173303e-5

bMin = -E/(kb*T)
bMax = 20.0

bMaxIndex = len(betas) -1

eq14 = []
eq14Neg = []
for b,beta in enumerate(betas):
    eq14.append(getEq14(beta,E,T,sab_water,A,alphas,len(betas),b))
    eq14Neg.append(getEq14(-beta,E,T,sab_water,A,alphas,len(betas),b))


invTotal = 1.0/(sum(eq14)+sum(eq14Neg))
eq14    = [ eq14Val * invTotal for eq14Val in eq14    ]
eq14Neg = [ eq14Val * invTotal for eq14Val in eq14Neg ]
plt.plot(betas,eq14,'ro-',markersize=2)
plt.plot([-b for b in betas],eq14Neg,'ro-',markersize=2)
plt.show()


