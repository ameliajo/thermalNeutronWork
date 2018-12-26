from plotSAB_help import *
import matplotlib.pyplot as plt
import sys
import subprocess
import matplotlib.colors as colors
import matplotlib.cm as cmx

##############################################################################
# READ IN DATA FROM NJOY-GENERATED FILE
##############################################################################
if len(sys.argv) > 1:
    if sys.argv[1] == 'njoy':
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/sabContin.txt','./'])
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/sabDelta.txt','./'])

def getLine(f):
    return [float(num) for num in f.readline().split()]

with open('sabContin.txt','r') as f:
    sabCONTIN = getLine(f)
    alphaVals = getLine(f)
    betaVals  = getLine(f)

with open('sabDelta.txt','r') as f:
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
ax.set_facecolor('xkcd:off white')
ax.set_facecolor('xkcd:light grey blue')
plt.show()




