import sys
sys.path.append('NJOY_LEAPR')
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as c
import matplotlib.cm as cmx
from generateNjoyInput import *
from getSAB import *
import numpy as np

# The point of this is program is to run the H in H2O LEAPR case, at T=296K,
# with its normal delta-function representation (phonon distribution is 
# approximated using delta functions at higher energy), as well as with a 
# similar triangle representation. LEAPR handles triangles different than delta
# functions, and so I'm looking at how changing the width of the triangles will
# impact the resultant S(a,b).


xs_bound = 2.0
kb = 8.61733e-5
T = 296.0
E = 1.0
mu_vec = list(np.linspace(-1,1,12))
Ep_vec = list(np.linspace(E,1.5,20))
Ep_vec = list(np.linspace(0.0,2.0,500))
A = 18.0


alphas = [1e-7+ 0.05*i for i in range(500)]
betas = [0.5*i for i in range(100)]

 
 # This is the part of the phonon distribution that is not usually approximated
 # using delta functions (lower E). 
continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]


oscE = [ 0.204,    0.4794   ]
oscW = [ 0.166667, 0.333333 ]

widths = list(range(2,12,2))
widths = [2] 

NJOY_LEAPR = False
NJOY_LEAPR = True
fullRedo = True
fullRedo = False
sabDELTA = getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    
sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) for width in widths]


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


for mu in mu_vec:
    sab_vec = []
    for Ep in Ep_vec:
        alpha = ( Ep + E - (2*mu*(Ep*E)**0.5) )/(A*kb*T)
        beta  = ( Ep - E )/(kb*T)

        a = findBounds(alphas,alpha)
        b = findBounds(betas,abs(beta))

        alphaL = alphas[a]
        alphaR = alphas[a+1]
        betaL = betas[b]
        betaR = betas[b+1]

        sab_betaL = interpolate(alphaL,alphaR,sabDELTA[a*len(betas)+b],sabDELTA[(a+1)*len(betas)+b],alpha)
        sab_betaR = interpolate(alphaL,alphaR,sabDELTA[a*len(betas)+b+1],sabDELTA[(a+1)*len(betas)+b+1],alpha)

        sabVal = interpolate(betaL,betaR,sab_betaL,sab_betaR,abs(beta))
        if (beta < 0):
            sabVal *= np.exp(abs(beta))

        sab_vec.append(sabVal)

    plt.plot(Ep_vec,sab_vec,label='mu = '+str(mu))
plt.legend(loc='best')
plt.show()






