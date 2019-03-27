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
    return None

def interpolate(x1,x2,y1,y2,x):
    m = (y2-y1)/(x2-x1)
    b = y2 - m*x2
    return m*x + b



def prepPlot(vec):
    cnorm = c.Normalize(vmin=0,vmax=len(vec)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = c.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vec))])
    colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,colorBarName):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel(colorBarName)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    plt.xlabel("E' (eV)")
    plt.ylabel("XS (b)")
    plt.show()

alphas = [1e-7+ 0.10*i for i in range(100)]
betas = [0.5*i for i in range(100)]
alphas = list(np.linspace(1e-7,80,300))
betas = list(np.linspace(0.0,50,800))
continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]


oscE = [ 0.204,    0.4794   ]
oscW = [ 0.166667, 0.333333 ]




def getXS_from_SAB(sab,alphas,betas,E,kb,T,Ep_vec,mu_vec,A):
    xs_each_mu = []
    for i,mu in enumerate(mu_vec):
        sab_vec = []
        xs_vec = []
        for Ep in Ep_vec:
            alpha = ( Ep + E - (2*mu*(Ep*E)**0.5) )/(A*kb*T)
            beta  = ( Ep - E )/(kb*T)

            a,b = findBounds(alphas,alpha), findBounds(betas,abs(beta))
            if a == None or b == None:
                xs_vec.append(0.0)
                continue
    
            alphaL, alphaR = alphas[a], alphas[a+1]
            betaL , betaR  = betas[b] , betas[b+1]

            sab_betaL = interpolate(alphaL,alphaR,sab[a*len(betas)+b],\
                                    sab[(a+1)*len(betas)+b],alpha)
            sab_betaR = interpolate(alphaL,alphaR,sab[a*len(betas)+b+1],\
                                     sab[(a+1)*len(betas)+b+1],alpha)

            sabVal = interpolate(betaL,betaR,sab_betaL,sab_betaR,abs(beta))
            #if(mu < 0.9999):
            #    print(Ep,mu,alpha,beta,a,b)
            #    print(sab_betaL,sabVal,sab_betaR)
            #    print()

            if (beta > 0): 
                sabVal *= np.exp(-beta)

            xs_vec.append((xs_bound/(2*kb*T))*(Ep/E)**0.5*sabVal)
            #xs_vec.append(sabVal)
        xs_each_mu.append(xs_vec)
    return xs_each_mu




if __name__=="__main__":
    plot_Ep_mu = False
    plot_Ep_Width = True
    plot_Ep_Width = False
    plot_Ep_Width_pct_errors = True
    plot_Ep_Width_pct_errors = False
    plot_Ep_errors = False
    plot_Ep_errors = True
    
    xs_bound = 20.449
    kb = 8.61733e-5
    T = 296.0
    E = 1.0
    A = 0.99917
    Ep_vec = list(np.linspace(0.0,1.5,1000))
    widths = list(range(2,12,2))

    NJOY_LEAPR = True
    fullRedo = True
    fullRedo = False
    sabDELTA   =  getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) \
                  for width in widths]


    if plot_Ep_mu:
        mu_vec = list(np.linspace(-1,1,64))
        #mu_vec = list(np.linspace(0.98,1,2))
        scalarMap, colorBar = prepPlot(mu_vec)
        xs_vec = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,mu_vec,A)
        for i in range(len(mu_vec)):
            plt.plot(Ep_vec,xs_vec[i],color=scalarMap.to_rgba(i))
            
        finishPlotting(colorBar,'mu')

    if plot_Ep_Width:
        mu = 0.9
        scalarMap, colorBar = prepPlot([0]+widths)

        xsDELTA = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0]
        xsCONTINS = [getXS_from_SAB(sabCONTIN,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0] \
                     for sabCONTIN in sabCONTINS]

        plt.plot(Ep_vec,xsDELTA,color=scalarMap.to_rgba(0),label='delta')
        for i,xsCONTIN in enumerate(xsCONTINS):
            plt.plot(Ep_vec,xsCONTIN,color=scalarMap.to_rgba(i+1),label='width = '+str(widths[i]))
 
        plt.legend(loc='best')
        plt.xlabel("E' (eV)")
        plt.ylabel("XS (b)")
        plt.show()

    if plot_Ep_Width_pct_errors:
        mu = 0.9
        scalarMap, colorBar = prepPlot([0]+widths)

        xsDELTA = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0]
        xsCONTINS = [getXS_from_SAB(sabCONTIN,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0] \
                     for sabCONTIN in sabCONTINS]

        for j,xsCONTIN in enumerate(xsCONTINS):
            #pctError = [100.0*(xsCONTIN[i]-xsDELTA[i])/xsDELTA[i] for i in range(len(xsDELTA))] \
            #        if abs(xsDELTA[i]) > 1e-6 else \
            #        [100.0*(xsCONTIN[i]-xsDELTA[i]) for i in range(len(xsDELTA))] 
            toPlot = []
            for i in range(len(xsDELTA)):
                if xsCONTIN[i] < 1e-12 and xsDELTA[i] < 1e-12: toPlot.append(0.0)
                else: toPlot.append(100.0*(xsCONTIN[i]-xsDELTA[i]) / \
                                   (0.5*(abs(xsCONTIN[i])+abs(xsDELTA[i]))))
                #else: toPlot.append(abs(xsCONTIN[i]-xsDELTA[i]))
                #toPlot.append((xsCONTIN[i]-xsDELTA[i]))

            plt.plot(Ep_vec,toPlot,color=scalarMap.to_rgba(j+1),label='width = '+str(widths[j]))
        plt.legend(loc='lower left')
        plt.xlabel("E' (eV)")
        plt.ylabel("Error (%)")
        plt.show()




    if plot_Ep_errors:
        mus = np.linspace(-1,1,21)
        scalarMap, colorBar = prepPlot([0]+mus)
        for j,mu in enumerate(mus):
            print(mu)
            xsDELTA = getXS_from_SAB(sabDELTA,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0]
            xsCONTINS = [getXS_from_SAB(sabCONTIN,alphas,betas,E,kb,T,Ep_vec,[mu],A)[0] \
                         for sabCONTIN in sabCONTINS]

            totalError = []
            for xsCONTIN in xsCONTINS:
                relError = [abs(xsCONTIN[i]-xsDELTA[i])/xsDELTA[i] for i in range(len(xsDELTA))]
                totalError.append(np.trapz(relError,x=Ep_vec))


            plt.plot(widths,totalError,color=scalarMap.to_rgba(j+1),label='mu = '+'%.3E'%mu)
            plt.plot(widths,totalError,color=scalarMap.to_rgba(j+1),marker='o')
        #plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel('triangle widths (# grid spaces)')
        plt.ylabel("Absolute Error (b)")
        ax = plt.gca()
        plt.colorbar(colorBar).ax.set_ylabel('mu')
        plt.show()






















