from plotSAB_help import *
import matplotlib.pyplot as plt
import sys
import subprocess
import matplotlib.colors as colors
import matplotlib.cm as cmx

ranNJOY = False
if len(sys.argv) > 1:
    if sys.argv[1] == 'njoy':
        ranNJOY = True


if ranNJOY:
    subprocess.run(['cp','/Users/amelia/NJOY2016/bin/sabContin.txt','./'])
    subprocess.run(['cp','/Users/amelia/NJOY2016/bin/sabDelta.txt','./'])


sabCONTIN = [float(num) for num in open('sabContin.txt','r').read().split()]
sabDELTA  = [float(num) for num in open('sabDelta.txt','r').read().split()]


A0 = 18.02
E = 0.01 
kbT = 0.025851
nalpha = 7
nbeta  = 11
alphaVals = [0.2,0.4,0.6,0.8,1.0,1.2,1.4]
betaVals  = [0,2,4,6,8,10,12,14,16,18,20]

alphaVals = [0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.8,1.0,1.2,1.4]
betaVals  = [0,2,4,6,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9.0,10,12,14,16,18,20]

betaVals  = [0,0.5,1.0,1.5,2,4,6,7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9.0,10,12,14,16,18,20]
alphaVals = [0.01,0.1,0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.8,1.0,1.2,1.4]

jet = cm = plt.get_cmap('hot') # hot, autumn, tab10
cnorm = colors.Normalize(vmin=0,vmax=len(alphaVals)+3)
scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=jet)


min,max = (-40,30)
step = 10
mymap = colors.LinearSegmentedColormap.from_list('funTestColors',[scalarMap.to_rgba(a) for a in range(len(alphaVals))])
#mymap = colors.LinearSegmentedColormap.from_list(scalarMap.to_rgba)
Z = [[0,0],[0,0]]
levels = range(min,max+step,step)
levels = alphaVals 
print(levels)
CS3 = plt.contourf(Z, levels, cmap=mymap)
plt.clf()


plotBetaForVariousAlpha(alphaVals,betaVals,sabDELTA,A0,E,kbT,scalarMap,'.',False)
plotBetaForVariousAlpha(alphaVals,betaVals,sabCONTIN,A0,E,kbT,scalarMap,'.',True)
#plotErrorBetaForVariousAlpha(alphaVals,betaVals,sabDELTA,sabCONTIN,A0,E,kbT,scalarMap)

ax = plt.gca()
cb = plt.colorbar(CS3) # using the colorbar info I got from contourf
cb.ax.set_ylabel('alpha values')


plt.title('S(a,b) values for water, generated with and without delta funcs')
ax.set_facecolor('xkcd:off white')
ax.set_facecolor('xkcd:light grey blue')

plt.show()




