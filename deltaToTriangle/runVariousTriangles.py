import sys
sys.path.append('NJOY_LEAPR')
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as c
import matplotlib.cm as cmx
from generateNjoyInput import *
from getSAB import *

# The point of this is program is to run the H in H2O LEAPR case, at T=296K,
# with its normal delta-function representation (phonon distribution is 
# approximated using delta functions at higher energy), as well as with a 
# similar triangle representation. LEAPR handles triangles different than delta
# functions, and so I'm looking at how changing the width of the triangles will
# impact the resultant S(a,b).


alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16, 18, 20]
 
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
sabDELTA = getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    
sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) for width in widths]

# These parameters are for making sure that we only consider logicat alpha, beta
# combinations. 
A0 = 18.02; E = 0.01; kbT = 0.025851

cMap = cmx.ScalarMappable(c.Normalize(0,10),plt.get_cmap('tab10')) #hot autumn tab10
colors = [cMap.to_rgba(i) for i in range(2*len(widths)+1)]

a = 9
plt_SAB_given_A(alphas,a,betas,sabDELTA,A0,E,kbT,colors[0],'.','with delta funcs',3)
for i in range(len(widths)):
    plt_SAB_given_A(alphas,a,betas,sabCONTINS[i],A0,E,kbT,colors[i+1],\
                    '.','triangle, width = '+ str(widths[i])+' spaces',3)



#plt.legend(loc='best')
ax = plt.gca()
LEAPR_type = "NJOY LEAPR" if NJOY_LEAPR else "MY LEAPR"
plt.title('S(a,b) values for water, delta vs. various triangles.\n Generated using '+LEAPR_type+', for alpha = '+str(alphas[a]))
ax.set_facecolor('xkcd:light grey blue')
ax.set_facecolor('xkcd:off white')
#plt.show()




NJOY_LEAPR = True
fullRedo = False
sabDELTA = getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,None,oscE,oscW)
    
sabCONTINS = [getSAB(alphas,betas,continRho,NJOY_LEAPR,fullRedo,width,oscE,oscW) for width in widths]

plt_SAB_given_A(alphas,a,betas,sabDELTA,A0,E,kbT,colors[0+len(widths)],'X','with delta funcs')
for i in range(len(widths)):
    plt_SAB_given_A(alphas,a,betas,sabCONTINS[i],A0,E,kbT,colors[i+1+len(widths)],\
                    'X','triangle, width = '+ str(widths[i])+' spaces',1)

plt.legend(loc='best')
plt.show()

