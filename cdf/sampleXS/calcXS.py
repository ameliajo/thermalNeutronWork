#from CDF_alpha import *
from CDF_beta import *
from random import random
import matplotlib.pyplot as plt
from rho import continRho
from getSAB import *
from alpha_CDF_no_bound_checking import *

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]



#for i in range(len(alphas_1)): assert(abs(alphas_1[i]-alphas_2[i]) < 1e-6)
#for i in range(len(betas_1)):  assert(abs(betas_1[i] -betas_2[i])  < 1e-6)
#for i in range(len(temps_1)):  assert(abs(temps_1[i] -temps_2[i])  < 1e-6)


def getAlphaMinMax(E,beta,T,A):
    kb = 8.6173303e-5
    return ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T ), \
           ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )

def getAlphaMinMaxIndices(E,beta,T,A,alphas):
    aMin,aMax = getAlphaMinMax(E,beta,T,A)
    aMin_index, aMax_index = 0, 0
    for a in range(len(alphas)):
        if alphas[a] > aMin and aMin_index == 0: aMin_index = a
        if alphas[a] <= aMax:                    aMax_index = a
    return aMin_index,aMax_index


alphas = alphas_2
betas = betas_2
temps = temps_2

t = 0

kb = 8.6173303e-5
T = temps[0]
#CDF_alpha = alphaCDF[t]
CDF_beta = betaCDF[t]
E = 1.0
A = 0.98

betaMin = -E/(kb*T)
betaMax = 20.0

N = 100

fullBetas = [-x for x in betas[::-1]] + betas[1:]
betas = [x for x in fullBetas if betaMin <=x <betaMax]


selectedAlphas = []
selectedBetas = []
for i in range(N):
    xsi_1 = random()
    index = 0
    beta = 0.0
    for i in range(len(CDF_beta)-1):
        if CDF_beta[i] <= xsi_1 < CDF_beta[i+1]:
            index = i
            fraction = (xsi_1 - CDF_beta[i]) / (CDF_beta[i+1] - CDF_beta[i])
            print(len(betas),len(CDF_beta))
            beta = (betas[i+1] - betas[i]) * fraction + betas[i]
            break
    

    sabsMINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=True,fullRedo=False,\
                width=None,oscE=oscE,oscW=oscW) for T in temps]
    CDF_alpha = PDF_CDF_at_various_temps(A,E,temps,beta,alphas,len(betas),index,sabsMINE)[t]


    aMin_index,aMax_index = getAlphaMinMaxIndices(E,beta,T,A,alphas)
    alphasToPlot = []
    H = []
    for a in range(aMin_index,aMax_index+1):
        alphasToPlot.append(alphas[a])
        H.append( (CDF_alpha[a]-CDF_alpha[aMin_index]) / \
                  (CDF_alpha[aMax_index]-CDF_alpha[aMin_index]) )



    xsi_2 = random()
    alpha = 0.0
    for i in range(len(CDF_alpha)-1):
        if CDF_alpha[i] <= xsi_2 < CDF_alpha[i+1]:
            fraction = (xsi_2 - CDF_alpha[i]) / (CDF_alpha[i+1] - CDF_alpha[i])
            alpha = (alphas[i+1] - alphas[i]) * fraction + alphas[i]
            break


    #print(alpha,beta)
    selectedAlphas.append(alpha)
    selectedBetas.append(beta)

plt.plot(selectedAlphas,selectedBetas,'ro')
plt.show()










