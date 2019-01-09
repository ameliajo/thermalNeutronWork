from getSAB import *
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot():
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    plt.title(title)
    ax.set_facecolor('xkcd:light grey blue') # off white
    plt.show()



def runAlphasBetas(opt0,opt1,alphas,betas,continRho,oscE,oscW,fullRedo=False):
    useNJOY1 = True if opt0[0] == 'njoy' else False if opt0[0] == 'mine' else None
    useNJOY2 = True if opt1[0] == 'njoy' else False if opt1[0] == 'mine' else None
    width0, width1 = opt0[1], opt1[1]

    assert useNJOY1 != None and useNJOY2 != None, \
           'Invalid identifiers. Use either "mine" or "njoy" for LEAPR type'
    assert useNJOY1 != useNJOY2 or width0 != width1,\
           'Comparing two identical cases'

    if useNJOY1 == useNJOY2:
        title  = 'NJOY' if useNJOY1 else 'my'
        title += ' LEAPR.\nUsing delta func. approx vs thin triangle approx.'
    else:
        assert width0 == width1, 'If comparing different widths, use same LEAPR'
        title  = 'my\nLEAPR vs NJOY LEAPR. Using '
        title += 'delta function approx.' if width0 == None else 'thin triangle approx.'
    A0 = 18.02; E = 0.01; kbT = 0.025851

    scalarMap, colorBar = prepPlot()
    ########################################################################
    # Generate specified S(a,b)
    SAB_1 = getSAB(alphas,betas,continRho,useNJOY1,fullRedo,width0,oscE,oscW)
    SAB_2 = getSAB(alphas,betas,continRho,useNJOY2,fullRedo,width1,oscE,oscW)
    ########################################################################
    # Compare S(a,b) values
    plotBetaForVariousAlpha(alphas,betas,SAB_1,A0,E,kbT,scalarMap,'.',False)
    plotBetaForVariousAlpha(alphas,betas,SAB_2,A0,E,kbT,scalarMap,'.',True)
    finishPlotting(colorBar,'S(a,b) for H in H2O, generated with '+title)
    ########################################################################
    # Compare S(a,b) errors
    plotErrorBetaForVariousAlpha(alphas,betas,SAB_1,SAB_2,A0,E,kbT,scalarMap)
    finishPlotting(colorBar,'S(a,b) error for H in H2O, generated with '+title)
    ########################################################################






if __name__=="__main__":
    alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 2.0, 3.0, 4.0 ]
    betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16,17, 18, 18.8, 18.9,19,19.1, 19.2, 20]
    continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
      .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
      .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
      .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
      .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
      .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
      .04244, .042, 0.0]

    oscE = [ 0.204,    0.4794   ] 
    oscW = [ 0.166667, 0.333333 ]

    fullRedo = True
    opt0,opt1 = ('mine',None),('mine',2) # Works
    #opt0,opt1 = ('njoy',None),('njoy',2) # Works
    #opt0,opt1 = ('mine',None),('njoy',None) # Works
    #opt0,opt1 = ('mine',2),('njoy',2) # Works
    #opt0,opt1 = ('mine',None),('njoy',2) # FAILS, diagonal comparison
    #opt0,opt1 = ('njoy',None),('mine',2) # FAILS, diagonal comparison
    #opt0,opt1 = ('mine',None),('mine',None) # FAILS, identical
    #opt0,opt1 = ('njoy',2),('njoy',2) # FAILS, identical
    #opt0,opt1 = ('test',None),('njoy',2) # FAILS, bad LEAPR identifier
    runAlphasBetas(opt0,opt1,alphas,betas,continRho,oscE,oscW,fullRedo)


