from plotSAB_help import *
import matplotlib.pyplot as plt
import sys
import subprocess
import matplotlib.colors as colors
import matplotlib.cm as cmx
from generateNjoyInput import *


##############################################################################
# READ IN DATA FROM NJOY-GENERATED FILE
##############################################################################
alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16, 18, 20]
continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]



if len(sys.argv) > 1:
    if sys.argv[1] == 'njoy':
        deltaName  = 'Delta'
        continName = 'Contin'
        # Generate S(a,b) using delta functions
        generateNjoyInput(deltaName,alphas,betas,continRho,True)
        runNJOY(deltaName)
        # Generate S(a,b) by approximating delta functions as a very thin triangle
        generateNjoyInput(continName,alphas,betas,getPhononDist(2,continRho),False)
        runNJOY(continName)



def getLine(f):
    return [float(num) for num in f.readline().split()]

with open('sabResults/sab_Contin.txt','r') as f:
    sabCONTIN = getLine(f)
    alphaVals = getLine(f)
    betaVals  = getLine(f)

with open('sabResults/sab_Delta.txt','r') as f:
    sabDELTA   = getLine(f)
    assert(getLine(f) == alphaVals)
    assert(getLine(f) == betaVals)



A0 = 18.02
E = 0.01 
kbT = 0.025851

cnorm = colors.Normalize(vmin=0,vmax=len(alphaVals)+3)
scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
        [scalarMap.to_rgba(a) for a in range(len(alphaVals))])
colorBar = plt.contourf([[0,0],[0,0]], alphaVals, cmap=mymap)
plt.clf()


plotBetaForVariousAlpha(alphaVals,betaVals,sabDELTA,A0,E,kbT,scalarMap,'.',False)
plotBetaForVariousAlpha(alphaVals,betaVals,sabCONTIN,A0,E,kbT,scalarMap,'.',True)
#plotErrorBetaForVariousAlpha(alphaVals,betaVals,sabDELTA,sabCONTIN,A0,E,kbT,scalarMap)


ax = plt.gca()
cb = plt.colorbar(colorBar) # using the colorbar info I got from contourf
cb.ax.set_ylabel('alpha values')
plt.title('S(a,b) values for water, generated with and without delta funcs')
#plt.title('Relative difference between S(a,b) generated with and without delta funcs')
ax.set_facecolor('xkcd:off white')
ax.set_facecolor('xkcd:light grey blue')
plt.show()




