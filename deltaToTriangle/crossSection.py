import sys
sys.path.append('NJOY_LEAPR')
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as c
import matplotlib.cm as cmx
from generateNjoyInput import *
from getSAB import *
import numpy as np

def findBounds(vec,val):
    for i in range(len(vec)-1):
        if vec[i] <= val <= vec[i+1]:
            return i
    print(vec,val)
    return None

def interpolate(x1,x2,y1,y2,x):
    m = (y2-y1)/(x2-x1)
    b = y2 - m*x2
    return m*x + b



def prepPlot(vec):
    cnorm = c.Normalize(vmin=0,vmax=len(vec)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    mymap = c.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vec))])
    colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,colorBarName):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel(colorBarName)
    ax.set_facecolor('xkcd:light grey blue') # off white
    plt.xlabel("E' (eV)")
    plt.ylabel("XS (b)")
    plt.show()

alphas = [1e-7+ 0.10*i for i in range(100)]
betas = [0.5*i for i in range(100)]
alphas = list(np.linspace(1e-7,15,100))
betas = list(np.linspace(0.0,50,300))
continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]


oscE = [ 0.204,    0.4794   ]
oscW = [ 0.166667, 0.333333 ]




def getXS_from_SAB(sab,alphas,betas,E,kb,T,Ep_vec,mu_vec):
    xs_each_mu = []
    for i,mu in enumerate(mu_vec):
        sab_vec = []
        xs_vec = []
        for Ep in Ep_vec:
            alpha = ( Ep + E - (2*mu*(Ep*E)**0.5) )/(A*kb*T)
            beta  = ( Ep - E )/(kb*T)

            a,b = findBounds(alphas,alpha), findBounds(betas,abs(beta))
    
            alphaL, alphaR = alphas[a], alphas[a+1]
            betaL , betaR  = betas[b] , betas[b+1]

            sab_betaL = interpolate(alphaL,alphaR,sab[a*len(betas)+b],\
                                    sab[(a+1)*len(betas)+b],alpha)
            sab_betaR = interpolate(alphaL,alphaR,sab[a*len(betas)+b+1],\
                                     sab[(a+1)*len(betas)+b+1],alpha)

            sabVal = interpolate(betaL,betaR,sab_betaL,sab_betaR,abs(beta))
            if (beta < 0): 
                sabVal *= np.exp(abs(beta))

            xs_vec.append((xs_bound/2*kb*T)*(Ep/E)**0.5*sabVal)
            #xs_vec.append(sabVal)
        xs_each_mu.append(xs_vec)
    return xs_each_mu




if __name__=="__main__":
    plot_Ep_mu = True
    plot_Ep_Width = False
    
    xs_bound = 20.449
    kb = 8.61733e-5
    T = 296.0
    E = 1.0
    A = 18.0
    Ep_vec = list(np.linspace(0.0,1.5,500))
    widths = list(range(2,12,2))

    NJOY_LEAPR = False
    fullRedo = False
    fullRedo = True
    sabDELTA   = getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) \
                  for width in widths]

    if plot_Ep_mu:
        mu_vec = list(np.linspace(-1,1,20))
        scalarMap, colorBar = prepPlot(mu_vec)
        xs_vec = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,mu_vec)
        for i in range(len(mu_vec)):
            plt.plot(Ep_vec,xs_vec[i],color=scalarMap.to_rgba(i))
        finishPlotting(colorBar,'mu')

    if plot_Ep_Width:
        mu = 0.5
        scalarMap, colorBar = prepPlot([0]+widths)

        xsDELTA = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,[mu])[0]
        xsCONTINS = [getXS_from_SAB(sabCONTIN,alphas,betas,E,kb,T,Ep_vec,[mu])[0] \
                     for sabCONTIN in sabCONTINS]

        plt.plot(Ep_vec,xsDELTA,color=scalarMap.to_rgba(0))
        for i,xsCONTIN in enumerate(xsCONTINS):
            plt.plot(Ep_vec,xsCONTIN,color=scalarMap.to_rgba(i+1))

        finishPlotting(colorBar,'triangleWidth')





















