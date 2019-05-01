import numpy as np
from math import factorial as f
import matplotlib.pyplot as plt 

alpha = 3
beta_i = 0.5
lambda_s = np.cosh(beta_i/2)/(np.sinh(beta_i/2)*beta_i)
tList = np.linspace(0,47,500)



def getExpansion(beta_i,tList,lambda_s,alpha):
    term1 = [ alpha*np.exp(-beta_i*0.5)*np.exp(-1j*beta_i*t) / \
              (2*beta_i*np.sinh(beta_i/2))+ \
              alpha*np.exp(beta_i*0.5)*np.exp(1j*beta_i*t) /\
              (-2*beta_i*np.sinh(-beta_i/2))\
              for t in tList]

    expansionVersion = [0.0]*len(term1)
    for t_i in range(len(tList)):
        for n in range(100):
            expansionVersion[t_i] += ((np.exp(-alpha*lambda_s)/f(n))*term1[t_i]**n)
    return expansionVersion

def getCalculated(beta_i,tList,alpha):
    gamma = [alpha/(2*beta_i*np.sinh(beta_i/2))  * (1-np.exp(-1j*beta_i*t))\
            *np.exp(-beta_i/2) + \
            alpha/(-2*beta_i*np.sinh(-beta_i/2)) * (1-np.exp(1j*beta_i*t))\
            *np.exp(beta_i/2) for t in tList]
    return [np.exp(-gamma[t]) for t in range(len(tList))]

def getRho(beta,delta,beta_i):
    return delta/(delta) * (-abs(beta-beta_i)+delta) if \
           (beta_i-delta<beta<beta_i+delta) else 0.0

def doTriangle(delta):
    betas = list(np.linspace(beta_i-delta,beta_i+delta,101))
    rho = [getRho(beta,delta,beta_i) for beta in betas]
    rhoArea = np.trapz(rho,x=betas)
    rho = [x/rhoArea for x in rho]
    
    P = [rho[b]/(2.0*betas[b]*np.sinh(betas[b]*0.5)) for b in range(len(betas))]
    integral = [alpha*np.trapz([P[b]*np.exp(-beta*0.5)*np.exp(-1j*beta*t) + \
                                P[b]*np.exp( beta*0.5)*np.exp( 1j*beta*t) \
                          for b,beta in enumerate(betas)],x=betas) for t in tList]


    triangleVersion = [0.0]*len(integral)
    for t_i in range(len(tList)):
        for n in range(50):
            triangleVersion[t_i] += ((np.exp(-alpha*lambda_s)/f(n))*integral[t_i]**n)
    return triangleVersion


plotActualValues = False
plotAbsErrors = True
plotAbsErrors = False 
plotPeakErrors = False


calculatedValue = getCalculated(beta_i,tList,alpha) if (plotAbsErrors or plotPeakErrors) else None
peakErrors = []

widths = [5e-2,1e-2,5e-3,1e-3,5e-4]
for width in widths:
    triangleExp = doTriangle(width)
    if plotActualValues:
        plt.plot(tList,triangleExp,label=str('%.0E'%width+' eV'))
    if plotAbsErrors:
        plt.plot(tList,[triangleExp[i]-calculatedValue[i] for i in range(len(tList))],label='%.1E'%width)
    if plotPeakErrors:
        peakErrors.append([triangleExp[i]-calculatedValue[i] for i in range(len(tList))])


if plotPeakErrors:
    chunk = int(len(tList)/4.5)
    bigDiffs_t = []
    bigDiffs_p = []
    for width_index,width in enumerate(widths):
        these_bigDiffs_t = []
        these_bigDiffs_p = []
        for i in range(4):
            current_t = tList[i*chunk:(i+1)*chunk]
            current_p = peakErrors[width_index][i*chunk:(i+1)*chunk]
            good_t, good_p = 0.0, 0.0
            for j,t in enumerate(current_t):
                if abs(current_p[j]) > abs(good_p):
                    good_p = current_p[j]
                    good_t = current_t[j]
            these_bigDiffs_t.append(good_t)
            these_bigDiffs_p.append(good_p)
            #plt.plot(tList[i*chunk:(i+1)*chunk],peakErrors[0][i*chunk:(i+1)*chunk])
        bigDiffs_t.append(these_bigDiffs_t)
        bigDiffs_p.append(these_bigDiffs_p)

plt.show()

#for i,width in enumerate(widths):
#    t = bigDiffs_t[i]
#    p = bigDiffs_p[i]
#    plt.plot(t,p,label=str('%.0E'%width+' eV'),marker='o')
#plt.xlabel('t [unitless]')
#plt.ylabel('exp [-gamma]')
#plt.show()

#for i in range(4):
#    t_s = []
#    p_s = []
#    for j,width in enumerate(widths):
#        t_s.append(bigDiffs_t[j][i])
#        p_s.append(bigDiffs_p[j][i])
#    plt.plot(widths,p_s,label=str('Location = '+'%.2f'%t_s[0]),marker='o')
#plt.xlabel('t [unitless]')
#plt.ylabel('exp [-gamma]')
#plt.show()


diff1 = [0.496,.0272,0.006,0.0002,0.0000723]
diff2 = [0.94,0.11,0.03,0.001,0.0002]
diff3 = [0.993,0.24,0.066,0.002,0.0006]
plt.plot(widths,diff1,marker='o',label='12 peak')
plt.plot(widths,diff2,marker='o',label='25 peak')
plt.plot(widths,diff3,marker='o',label='37 peak')
plt.xlabel('Triangle Width [unitless]')
plt.ylabel('Absolute difference in exp[-gamma]')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
plt.show()


#expansion = getExpansion(beta_i,tList,lambda_s,alpha)
#plt.plot(tList,expansion,label="expansion")
#calculated = getCalculated(beta_i,tList,alpha)
#plt.plot(tList,calculated,label="calculated")

#plt.ylabel('exp [-gamma]')
#plt.xlabel('t [unitless]')
#if plotActualValues: plt.legend(loc='upper right')
#if plotAbsErrors:    plt.legend(loc='lower left')
#plt.show()

