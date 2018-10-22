import matplotlib.pyplot as plt
from results_test_09_300K import *
from math import exp




def getEq14(beta,E,T,Sab,A,alphaVec,lenBeta,b):
    kb = 8.6173303e-5
    alpha_min = ( (E)**0.5 - (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    alpha_max = ( (E)**0.5 + (E+(beta)*kb*T)**0.5 )**2 / ( A*kb*T )
    validAlpha = [ (i,alphaVec[i]) for i in range(len(alphaVec)) if \
                   (alphaVec[i] >= alpha_min and alphaVec[i] <= alpha_max) ]

    numerator, denominator = 0.0, 0.0
    asymConverter = exp(-(beta)/2)
    for i_alpha in validAlpha:
        i = i_alpha[0] 
        numerator += Sab[i*lenBeta+b]
    if beta < 0.0 and round(beta) == -19: print(asymConverter,numerator)
    if beta < 0.0 and round(beta) == -20: print()
    if beta < 0.0 and round(beta) == -20: print(asymConverter,numerator)
    if beta < 0.0 and round(beta) == -22: print()
    if beta < 0.0 and round(beta) == -22: print(asymConverter,numerator)
    numerator *= asymConverter

    return numerator 



A = 18.02   # H2O = 1.01 + 1.01 + 16.00 = 18.02 
E = 1       # 1 eVc
T = 300     # 300 K
kb = 8.6173303e-5

bMin = -E/(kb*T)
bMax = 20.0


eq14 = []
eq14_earlier = []
for b,beta in enumerate(betas):
    eq14.append(getEq14(beta,E,T,sab_water,A,alphas,len(betas),b))
    eq14_earlier.append(getEq14(-beta,E,T,sab_water,A,alphas,len(betas),b))


invTotal = 1.0/(sum(eq14)+sum(eq14_earlier))
eq14 = [ eq14Val * invTotal for eq14Val in eq14 ]
eq14_earlier = [ eq14Val * invTotal for eq14Val in eq14_earlier ]

plt.plot([-b for b in betas],eq14_earlier,'ro-',markersize=2)
plt.plot(betas,eq14,'ro-',markersize=2)
plt.show()
"""
"""






