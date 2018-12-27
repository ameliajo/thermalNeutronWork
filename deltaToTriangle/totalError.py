from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as c
import matplotlib.cm as cmx
from generateNjoyInput import *
import sys

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




# Define all the different triangle widths we're going to be looking at
widths = list(range(2,12,2))
deltaName = 'deltaInput'
fileNames = ['triangleOfWidth'+str(width) for width in widths]

# If necessary, we generate NJOY input files and give them to NJOY. 
if len(sys.argv) > 1:
    if sys.argv[1] == 'njoy':
        generateNjoyInput(deltaName,alphas,betas,continRho,True)
        runNJOY(deltaName)

        for i,width in enumerate(widths):
            generateNjoyInput(fileNames[i],alphas,betas,\
                              getPhononDist(width,continRho),False)
            runNJOY(fileNames[i])


def getLine(f):
    return [float(num) for num in f.readline().split()]


# Read in the S(a,b) output that was generated from NJOY. Note that currently
# the S(a,b) data lives in the sabResults/ subdirectory. 
sabCONTINS = []
for fileName in fileNames:
    with open('sabResults/sab_'+str(fileName)+'.txt','r') as f:
        sabCONTINS.append(getLine(f))
        assert(alphas == getLine(f)) # Sanity check, hoping we 
        assert(betas  == getLine(f)) # generated our inputs okay

with open('sabResults/sab_'+str(deltaName)+'.txt','r') as f:
    sabDELTA = getLine(f)
    assert(alphas == getLine(f))
    assert(betas  == getLine(f))



#cMap = cmx.ScalarMappable(c.Normalize(0,10),plt.get_cmap('tab10')) #hot autumn tab10
#colors = [cMap.to_rgba(i) for i in range(len(widths)+1)]

error = []
for i in range(len(widths)):
    c = sabCONTINS[i]
    d = sabDELTA
    error.append(sum([abs(c[i]-d[i]) for i in range(len(c))]))
    print(i,error[-1])


plt.plot(widths,error,marker='o')
plt.xlabel('# of spaces spanned by triangle')
plt.ylabel('Total absolute error')
plt.title('Total absolute error between S(a,b) generated using delta\nfunctions vs. S(a,b) generated using triangles of various widths')
plt.show()

#plt.legend(loc='best')
#ax = plt.gca()
#plt.title('S(a,b) values for water, generated with delta,\n and with triangles of various widths, for alpha = '+str(alphas[a]))
#ax.set_facecolor('xkcd:light grey blue')
#ax.set_facecolor('xkcd:off white')


#plt.show()

