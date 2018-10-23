from math import pi
from math import exp
import matplotlib.pyplot as plt
import numpy as np


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )


A = 18.02   # H2O = 1.01 + 1.01 + 16.00 = 18.02 
A = 1.0   # H2O = 1.01 + 1.01 + 16.00 = 18.02 
E = 1.0     # 1 eVc
T = 300     # 300 K
kb = 8.6173303e-5


bMin = -E/(kb*T)
bMax = 20.0


beta = np.linspace(-38,bMax,59)
# make sure that beta grid has 0 in it
numAlpha = 100.0
eq14 = []
asym = []
for b in beta:
    
    alpha_min = ( (E)**0.5 - (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
    alpha_max = ( (E)**0.5 + (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
    

    alpha = np.linspace(alpha_min,alpha_max,numAlpha)
    eq14.append(sum([calcAsym(a,b)*(alpha_max-alpha_min)/numAlpha for a in alpha]))

plt.plot(beta,eq14)
plt.show()

"""
#invIntegral = 1.0/np.trapz(eq14, x=beta)
#eq14 = [x * invIntegral for x in eq14]
b = -38.0
alpha_min = ( (E)**0.5 - (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
alpha_max = ( (E)**0.5 + (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
alpha = np.linspace(alpha_min,alpha_max,10)
for a in alpha:
    asym.append(calcAsym(a,b)*(alpha_max-alpha_min)/10.0)
for a in asym:
    print(a)
print()
print(sum(asym))
plt.plot(alpha,asym)
plt.show()

"""
    
    
    
#width = (alpha_max-alpha_min)/1000.0

#invTotal = 1.0/(sum(asym)*width)
#invTotal = 1.0/(sum(eq14)*width)
#eq14 = [entry * invTotal for entry in eq14]
#asym = [entry * invTotal for entry in asym]








