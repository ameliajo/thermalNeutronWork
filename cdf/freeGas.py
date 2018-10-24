from math import pi
from math import exp
import matplotlib.pyplot as plt
import numpy as np
from math import ceil


def calcAsym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return (1.0/(4.0*pi*alpha)**0.5) * exp( - (alpha+beta)**2 / (4.0*alpha) )

def calcSym(alpha,beta):
    if alpha == 0: alpha = 0.0001
    return exp(beta*0.5)*calcAsym(alpha,beta)


A = 1.0     # H2O = 1.01 + 1.01 + 16.00 = 18.02 
E = 1.0     # 1 eVc
T = 1000     # 300 K
kb = 8.6173303e-5

bMin = -E/(kb*T)
bMax = 20.0

beta = np.linspace(ceil(bMin),bMax,bMax-bMin+1)
# make sure that beta grid has 0 in it

bFile = open("betaFreeGas.py", "w")
#bFile.write("betas = "+str([float("%.5f"%round(b,6)) for b in beta]))
bFile.write("betas = "+str([b for b in beta]))
bFile.close()

aFile = open("alphaFreeGas.py", "w")
aFile.write("alphas = [")

sabFile = open("sabFreeGas.py", "w")
sabFile.write("sab = [")

sabSymFile = open("sabFreeGasSym.py", "w")
sabSymFile.write("sabSym = [")



numAlpha = 100.0
eq14 = []
for b in beta:
    
    alpha_min = ( (E)**0.5 - (E+b*kb*T)**0.5 )**2 / ( A*kb*T )
    alpha_max = ( (E)**0.5 + (E+b*kb*T)**0.5 )**2 / ( A*kb*T )

    alpha = np.linspace(alpha_min,alpha_max,numAlpha)
    if b == beta[-1]:
        #aFile.write(str([float("%.5f"%round(a,6)) for a in alpha])+"]")
        aFile.write(str([a for a in alpha])+"]")
    else:
        #aFile.write(str([float("%.5f"%round(a,6)) for a in alpha])+", ")
        aFile.write(str([a for a in alpha])+", ")
    sabVec = [calcAsym(a,b) for a in alpha]
    sabSymVec = [calcSym(a,b) for a in alpha]

    if b == beta[-1]:
        #sabFile.write(str([float("%.5f"%round(sab,6)) for sab in sabVec])+"]")
        sabFile.write(str([sab for sab in sabVec])+"]")
    else:
        #sabFile.write(str([float("%.5f"%round(sab,6)) for sab in sabVec])+", ")
        sabFile.write(str([sab for sab in sabVec])+", ")

    if b == beta[-1]: sabSymFile.write(str([sab for sab in sabSymVec])+"]")
    else: sabSymFile.write(str([sab for sab in sabSymVec])+", ")


    #eq14.append(sum([calcAsym(a,b)*(alpha_max-alpha_min)/numAlpha for a in alpha]))
    eq14.append(np.trapz(sabVec,x=alpha))
    #eq14.append(exp(-b/2)*(np.trapz(sabSymVec,x=alpha)))
    if ( b <= 0 ) : 
        print(eq14[-1])
    
    
aFile.close()

invIntegral = 1.0/np.trapz(eq14, x=beta)
eq14 = [ x * invIntegral for x in eq14 ]

f = plt.figure(1)
plt.plot(beta,eq14)
f.show()
hold = input()
"""


"""


