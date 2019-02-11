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


alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 7.95,8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16, 18, 20]
betas = [7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 7.95,8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0]
betas = [7.0 + 0.01*i for i in range(200)]
betas = [7.5 + 0.02*i for i in range(50)]
#betas = [7.0 + i for i in range(3)]
#alphas = [1.0*i for i in range(5)]

#alphas = [0.5]

 
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

NJOY_LEAPR = False
fullRedo = False
fullRedo = True
sabDELTA = getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    
sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) for width in widths]

# Sum over all values in the sab list, take difference
error = [ sum( [abs(sabCONTIN[j]-sabDELTA[j]) for j in range(len(sabDELTA))] ) \
          for sabCONTIN in sabCONTINS ]



def integrateAcrossBeta(sabCONTIN,sabDELTA,a,betas):
    val = 0.0
    for b in range(len(betas)-1):
        val_left  = abs(sabCONTIN[a*len(betas)+b]-sabDELTA[a*len(betas)+b])
        val_right = abs(sabCONTIN[a*len(betas)+b+1]-sabDELTA[a*len(betas)+b+1])
        val += (val_left+val_right)*0.5*(betas[b+1]-betas[b])
    return val

error = []
for sabCONTIN in sabCONTINS:
    errorVal = 0.0
    for a in range(len(alphas)-1):
       val_left = integrateAcrossBeta(sabCONTIN,sabDELTA,a,betas)
       val_right = integrateAcrossBeta(sabCONTIN,sabDELTA,a+1,betas)
       errorVal += (val_left + val_right)*0.5*(alphas[a+1]-alphas[a])
    error.append(errorVal)


widths = [0.00255*w for w in widths]
plt.xticks(np.arange(min(widths), max(widths)+0.00255,2*0.00255))
plt.plot(widths,error,marker='o')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel('Triangle Width (eV)')
plt.ylabel('Total Difference')
#plt.title('Total absolute error between S(a,b) generated using delta\nfunctions vs. S(a,b) generated using triangles of various widths')
plt.show()



