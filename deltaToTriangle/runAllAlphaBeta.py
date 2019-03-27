from getSAB import *
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot():
    cnorm = colors.Normalize(vmin=0,vmax=2*len(alphas)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    #plt.title(title)
    #ax.set_facecolor('xkcd:light grey blue') # off white
    plt.yscale('log')
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
    A0 = 1.0; E = 0.50; kbT = 0.025851

    scalarMap, colorBar = prepPlot()
    ########################################################################
    # Generate specified S_nonsym(a,-b)
    SAB_1 = getSAB(alphas,betas,continRho,useNJOY1,fullRedo,width0,oscE,oscW)
    SAB_2 = getSAB(alphas,betas,continRho,useNJOY2,fullRedo,width1,oscE,oscW)
    #print(betas_negative)

    for a in range(len(alphas)):
        for b in range(len(betas)):
            SAB_1[a*len(betas)+b] *= np.exp(-betas[b]/2)
            SAB_2[a*len(betas)+b] *= np.exp(-betas[b]/2)






    ########################################################################
    # Compare S(a,b) values
    plotBetaForVariousAlpha(alphas,betas,SAB_1,A0,E,kbT,scalarMap,'.',True)
    plotBetaForVariousAlpha(alphas,betas,SAB_2,A0,E,kbT,scalarMap,'.',False)
    finishPlotting(colorBar,'S(a,b) for H in H2O, generated with '+title)
    ########################################################################
    # Compare S(a,b) errors (percent)
    #plotErrorBetaForVariousAlpha(alphas,betas,SAB_1,SAB_2,A0,E,kbT,scalarMap)
    #finishPlotting(colorBar,'S(a,b) error for H in H2O, generated with '+title)
    ########################################################################
    # Compare S(a,b) errors (absolute)
    #plotAbsoluteErrorBetaForVariousAlpha(alphas,betas,SAB_1,SAB_2,A0,E,kbT,scalarMap)
    #finishPlotting(colorBar,'S(a,b) error for H in H2O, generated with '+title)
    ########################################################################








if __name__=="__main__":
    alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 2.0, 3.0, 4.0, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    alphas = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9,1.0, 1.1, 1.2,1.3,1.4,1.5,1.6,1.7]
    betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.3, 7.4, 7.44, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16,17, 18, 18.7, 18.75, 18.8, 18.85,18.9,18.95,19,19.1, 19.2, 20,22,24,26,26.6,26.7,26.8,28,30,34,36,37.4,37.5,37.6,38]

    alphas = [.01008, .015, .0252, .033, 0.050406, .0756, 0.100812, 0.151218, 0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 0.504060, 0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 0.976435, 1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 1.873790, 2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 3.578320, 3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 6.816500, 7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 12.96740, 14.21450, 15.58150, 17.07960, 18.72080] 

    betas = [0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 0.065750, .0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 1.048440, 1.129090, 1.209740, 1.290390, 1.371040, 1.451690, 1.532340, 1.612990, 1.693640, 1.774290, 1.854940, 1.935590, 2.016240, 2.096890, 2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 2.580790, 2.669500, 2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 3.595360, 3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 9.637220, 10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 16, 16.5, 17, 17.5, 18, 18.5, 18.6, 18.7, 18.8, 18.9, 19, 19.1, 19.2, 19.3, 19.4, 19.5, 20, 20.5 ]

    betas = [0.03*i for i in range(1000)]
    alphas= [0.01*i for i in range(1,150)]

    alphas= [0.1*i for i in range(1,150)]
    betas = [0.1001*i for i in range(300)]
    betas = list(np.linspace(0,30,100))
    betas += [7.8+0.01*i for i in range(50)]
    betas += [15.8+0.01*i for i in range(50)]
    betas += [18.7+0.01*i for i in range(200)]
    betas += [26.5+0.01*i for i in range(50)]
    #betas += list(np.linspace(18.701,18.9,40))
    alphas = np.linspace(0.001,30,50)
    #alphas = [10.0,200]

    betas.sort()



    continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
      .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
      .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
      .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
      .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
      .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
      .04244, .042, 0.0]

    oscE = [ 0.204,    0.4794   ] 
    oscW = [ 0.166667, 0.333333 ]
    fullRedo = False
    fullRedo = True
    #opt0,opt1 = ('mine',None),('mine',2) # Works
    opt0,opt1 = ('njoy',None),('njoy',2) # Works
    #opt0,opt1 = ('njoy',None),('njoy',6) # Works
    #opt0,opt1 = ('njoy',None),('mine',None) # Works
    #opt0,opt1 = ('mine',None),('njoy',None) # Works
    #opt0,opt1 = ('mine',2),('njoy',2) # Works
    #opt0,opt1 = ('mine',None),('njoy',2) # FAILS, diagonal comparison
    #opt0,opt1 = ('njoy',None),('mine',2) # FAILS, diagonal comparison
    #opt0,opt1 = ('mine',None),('mine',None) # FAILS, identical
    #opt0,opt1 = ('njoy',2),('njoy',2) # FAILS, identical
    #opt0,opt1 = ('test',None),('njoy',2) # FAILS, bad LEAPR identifier
    runAlphasBetas(opt0,opt1,alphas,betas,continRho,oscE,oscW,fullRedo)


