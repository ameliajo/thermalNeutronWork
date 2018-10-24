from math import pi
from math import exp
import matplotlib.pyplot as plt
import numpy as np
from math import ceil


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )


A = 1.0     # H2O = 1.01 + 1.01 + 16.00 = 18.02 
E = 1.0     # 1 eVc
T = 1000     # 300 K
kb = 8.6173303e-5


bMin = -E/(kb*T)
bMax = 20.0
print(bMin,bMax)

beta = np.linspace(-38,bMax,59)
beta = np.linspace(ceil(bMin),bMax,bMax-bMin+1)
# make sure that beta grid has 0 in it
numAlpha = 100.0
eq14 = []
asym = []
for b in beta:
    
    alpha_min = ( (E)**0.5 - (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
    alpha_max = ( (E)**0.5 + (E+b*kb*T)**0.5 )**2 / ( A*kb*T )

    alpha = np.linspace(alpha_min,alpha_max,numAlpha)
    #eq14.append(sum([calcAsym(a,b)*(alpha_max-alpha_min)/numAlpha for a in alpha]))
    sabVec = [calcAsym(a,b) for a in alpha]
    eq14.append(np.trapz(sabVec,x=alpha))
    

invIntegral = 1.0/np.trapz(eq14, x=beta)
eq14 = [ x * invIntegral for x in eq14 ]

f = plt.figure(1)
plt.plot(beta,eq14)
f.show()
hold = input()




