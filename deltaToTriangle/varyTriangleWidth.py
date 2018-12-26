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
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/cw2.txt','./'])
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/cw4.txt','./'])
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/cw6.txt','./'])
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/cw8.txt','./'])
        subprocess.run(['cp','/Users/amelia/NJOY2016/bin/sabDelta.txt','./'])

def getLine(f):
    return [float(num) for num in f.readline().split()]

with open('cw2.txt','r') as f:
    sabCONTIN2 = getLine(f)
    alphaVals = getLine(f)
    betaVals  = getLine(f)

with open('cw4.txt','r') as f:
    sabCONTIN4 = getLine(f)
    assert(getLine(f) == alphaVals)
    assert(getLine(f) == betaVals)

with open('cw6.txt','r') as f:
    sabCONTIN6 = getLine(f)
    assert(getLine(f) == alphaVals)
    assert(getLine(f) == betaVals)

with open('cw8.txt','r') as f:
    sabCONTIN8 = getLine(f)
    assert(getLine(f) == alphaVals)
    assert(getLine(f) == betaVals)

with open('sabDelta.txt','r') as f:
    sabDELTA   = getLine(f)
    assert(getLine(f) == alphaVals)
    assert(getLine(f) == betaVals)



A0 = 18.02
E = 0.01 
kbT = 0.025851

cnorm = colors.Normalize(vmin=0,vmax=6)
scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10

alphaIndex = 5

plotBetaForGivenAlpha(alphaVals,alphaIndex,betaVals,sabDELTA,A0,E,kbT,scalarMap.to_rgba(0),'.','with delta funcs')
plotBetaForGivenAlpha(alphaVals,alphaIndex,betaVals,sabCONTIN2,A0,E,kbT,scalarMap.to_rgba(1),'.','triangle, width = 2 spaces')
plotBetaForGivenAlpha(alphaVals,alphaIndex,betaVals,sabCONTIN4,A0,E,kbT,scalarMap.to_rgba(2),'.','triangle, width = 4 spaces')
plotBetaForGivenAlpha(alphaVals,alphaIndex,betaVals,sabCONTIN6,A0,E,kbT,scalarMap.to_rgba(3),'.','triangle, width = 6 spaces')
plotBetaForGivenAlpha(alphaVals,alphaIndex,betaVals,sabCONTIN8,A0,E,kbT,scalarMap.to_rgba(4),'.','triangle, width = 8 spaces')


plt.legend(loc='best')
ax = plt.gca()
plt.title('S(a,b) values for water, generated with delta,\n and with triangles of various widths, for alpha = '+str(alphaVals[alphaIndex]))
ax.set_facecolor('xkcd:light grey blue')
ax.set_facecolor('xkcd:off white')
plt.show()




