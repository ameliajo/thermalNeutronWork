import numpy as np
import matplotlib.pyplot as plt
from math import pi
from getSAB import *


import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20c')) #hot autumn tab10
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    #plt.title(title)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    #ax.set_facecolor('xkcd:very light blue') # off white
    #plt.yscale('log')
    plt.show()




continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]

alphas = list(np.linspace(0.01,60,500))
betas = list(np.linspace(0.0,40,100))

n_alpha = len(alphas)
n_beta = len(betas)


oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]

fullRedo = False
fullRedo = True
runNJOY = False
width = None 
temps = [296.0,475.0,650.0,825.0,1000.0]   # Water
temps = list(np.linspace(296,1000,20))


def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]

def calcAsym(alpha,beta):
    return (1.0/(4.0*pi*alpha)**0.5) * np.exp( - (alpha+beta)**2 / (4.0*alpha) )

def calcSym(alpha,beta):
    return np.exp(beta*0.5)*calcAsym(alpha,beta)

def getAlphaMinMax(E,beta,kb,T,A):
    aMin = ( (E)**0.5 - (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    aMax = ( (E)**0.5 + (E+beta*kb*T)**0.5 )**2 / ( A*kb*T )
    return aMin,aMax

def getAlphaRange(E,beta,T,A,alphas):
    kb = 8.6173303e-5
    aMin,aMax = getAlphaMinMax(E,beta,kb,T,A)
    alphas = [x for x in alphas if aMin <= x <= aMax]
    return alphas

def calcAlphaPDF(A,E,t,beta,alphas,n_beta,index,sabs):
    denominator = 0.0
    for a in range(len(alphas)-1):
        sabL = getSABval(sabs[t],a,index,n_beta)
        sabR = getSABval(sabs[t],a+1,index,n_beta)
        denominator += (sabL+sabR)*0.5*(alphas[a+1]-alphas[a])

    eq15 = [0.0]*len(alphas)
    for a in range(len(eq15)):
        eq15[a] = getSABval(sabs[t],a,index,n_beta)/denominator
    return eq15

def calcAlphaCDF(eq15,alphas):
    eq16 = [0.0]*len(eq15)
    for a in range(1,len(alphas)):
        eq16[a] = eq16[a-1]+(eq15[a]+eq15[a-1])*0.5*(alphas[a]-alphas[a-1])
    return eq16

def PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabs):
    alpha_vecs, eq15_vecs, eq16_vecs = [], [], []
    for t in range(len(temps)):
        alphas = getAlphaRange(E,beta,temps[t],A,alphas)
        eq15 = calcAlphaPDF(A=A,E=E,t=t,beta=beta,alphas=alphas,n_beta=n_beta,index=index,sabs=sabs)
        alpha_vecs.append(alphas)
        eq15_vecs.append(eq15)
        eq16 = calcAlphaCDF(eq15,alphas)
        eq16_vecs.append(eq16)

    return alpha_vecs,eq15_vecs,eq16_vecs



A = 0.98
E = 1.0 

# A with 18 and E with 10 looks great
beta = 10
index = 63
index = 50
index = 25
beta = betas[index]
print("BETA",beta)
fullRedo = True
fullRedo = False

#sabsMINE = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=runNJOY,fullRedo=False,\
#            width=width,oscE=oscE,oscW=oscW) for T in temps]
sabsNJOY = [getSAB(alphas,betas,T,continRho,NJOY_LEAPR=True,fullRedo=fullRedo,\
            width=width,oscE=oscE,oscW=oscW) for T in temps]



alpha_vecs, alphaPDFs, alphaCDFs = PDF_CDF_at_various_temps(A,E,temps,beta,alphas,n_beta,index,sabsNJOY)
b_index = index

alphaValsToAimFor = list(np.linspace(5,35,20))
scalarMap, colorBar = prepPlot(alphaValsToAimFor)
for a,alphaValToAimFor in enumerate(alphaValsToAimFor):
    cdfVals = []
    for t in range(len(temps)):
        aPrime = 0
        for i in range(len(alpha_vecs[t])-1):
            if alpha_vecs[t][i] <= alphaValToAimFor <= alpha_vecs[t][i+1]: aPrime = i
        #pdfVals.append(alphaPDFs[t][aPrime])
        cdfVals.append(alphaCDFs[t][aPrime])
    plt.plot(temps,cdfVals,color=scalarMap.to_rgba(a))
    plt.plot(temps,cdfVals,'x',markersize=2,color=scalarMap.to_rgba(a))

plt.xlabel("Temperature (K)")
plt.ylabel("CDF")
finishPlotting(colorBar,"")























