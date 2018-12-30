import sys
sys.path.append('../')
import subprocess
from makeTest09Rho import *
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


# Things to plot: 
toPlot = ['delta','2']


def getLine(f):
    return [float(num) for num in f.readline().split()]


alphaVals = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
betaVals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16, 18, 20]

oscE = [ 0.204,    0.4794   ] 
oscW = [ 0.166667, 0.333333 ]



if len(sys.argv) > 1:
    if sys.argv[1] == 'partial' or sys.argv[1] == 'full':
        deltaName = 'sab_Delta.txt'
        writeRho('inputVals.txt',continRho,alphaVals,betaVals,oscE,oscW)
        if sys.argv[1] == 'full':
            subprocess.run(['g++','-std=c++14','deltaFuncLEAPR.cpp'])
        subprocess.run(['./a.out'])
        subprocess.run(['mv','inputVals.txt','alphaBetaRhoInputs/'+deltaName])
        subprocess.run(['mv','outputSAB.txt','sabResults/'+deltaName])

        continName = 'sab_Contin.txt'
        writeRho('inputVals.txt',getPhononDist(2,continRho),alphaVals,betaVals,[],[])
        if sys.argv[1] == 'full':
            subprocess.run(['g++','-std=c++14','deltaFuncLEAPR.cpp'])
        subprocess.run(['./a.out'])
        subprocess.run(['mv','inputVals.txt','alphaBetaRhoInputs/'+continName])
        subprocess.run(['mv','outputSAB.txt','sabResults/'+continName])


with open('sabResults/sab_Delta.txt','r') as f:
    sabDELTA  = getLine(f)

with open('sabResults/sab_Contin.txt','r') as f:
    sabCONTIN  = getLine(f)



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




#subprocess.run(['rm','a.out'])
subprocess.run(['rm','-rf','__pycache__'])



