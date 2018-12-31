
from getSAB import *
from plotSAB_help import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def finishPlotting(colorBar,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    plt.title(title)
    ax.set_facecolor('xkcd:light grey blue') # off white
    plt.show()




def runAlphasBetas(option):
    # Option
    # (1)  My  LEAPR delta function vs NJOY LEAPR delta function
    # (2)  My  LEAPR thin triangle  vs NJOY LEAPR thin triangle
    # (3)  My  LEAPR delta function vs  MY  LEAPR thin triangle
    # (4) NJOY LEAPR delta function vs NJOY LEAPR thin triangle

    alphas = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
    betas = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 9, 10, 12, 14, 16, 18, 20]
    continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
      .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
      .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
      .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
      .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
      .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
      .04244, .042, 0.0]

    oscE = [ 0.204,    0.4794   ] 
    oscW = [ 0.166667, 0.333333 ]


    A0 = 18.02
    E = 0.01 
    kbT = 0.025851

    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()

    fullRedo = False

    # (1)  My  LEAPR delta function vs NJOY LEAPR delta function
    # (2)  My  LEAPR thin triangle  vs NJOY LEAPR thin triangle
    # (3)  My  LEAPR delta function vs  MY  LEAPR thin triangle
    # (4) NJOY LEAPR delta function vs NJOY LEAPR thin triangle

    if option == 1:
        ########################################################################
        # COMPARE MY LEAPR DELTA FUNCS AGAINST NJOY LEAPR DELTA FUNCS
        ########################################################################
        MY_DELTA_SAB   = getSAB(alphas,betas,continRho,NJOY_LEAPR=False,\
                                fullRedo=fullRedo,width=None,oscE=oscE,oscW=oscW)
        NJOY_DELTA_SAB = getSAB(alphas,betas,continRho,NJOY_LEAPR=True,\
                                fullRedo=fullRedo,width=None,oscE=oscE,oscW=oscW)
        ########################################################################
        # Compare S(a,b) values
        plotBetaForVariousAlpha(alphas,betas,MY_DELTA_SAB,A0,E,kbT,scalarMap,'.',False)
        plotBetaForVariousAlpha(alphas,betas,NJOY_DELTA_SAB,A0,E,kbT,scalarMap,'.',True)
        finishPlotting(colorBar,'S(a,b) for H in H2O, generated using my LEAPR vs. NJOY LEAPR\nUsing delta function approx.')
        ########################################################################
        # Compare S(a,b) errors
        plotErrorBetaForVariousAlpha(alphas,betas,MY_DELTA_SAB,NJOY_DELTA_SAB,A0,E,kbT,scalarMap)
        finishPlotting(colorBar,'S(a,b) error for H in H2O, generated using my LEAPR vs. NJOY LEAPR\nUsing delta function approx.')
        ########################################################################


    if option == 2:
        ########################################################################
        # COMPARE MY LEAPR CONTIN2 FUNCS AGAINST NJOY LEAPR CONTIN2 FUNCS
        ########################################################################
        MY_CONTIN2_SAB   = getSAB(alphas,betas,continRho,NJOY_LEAPR=False,\
                                fullRedo=fullRedo,width=2,oscE=oscE,oscW=oscW)
        NJOY_CONTIN2_SAB = getSAB(alphas,betas,continRho,NJOY_LEAPR=True,\
                                fullRedo=fullRedo,width=2,oscE=oscE,oscW=oscW)
        ########################################################################
        # Compar S(a,b) values
        plotBetaForVariousAlpha(alphas,betas,MY_CONTIN2_SAB,A0,E,kbT,scalarMap,'.',False)
        plotBetaForVariousAlpha(alphas,betas,NJOY_CONTIN2_SAB,A0,E,kbT,scalarMap,'.',True)
        finishPlotting(colorBar,'S(a,b) for H in H2O, generated using my LEAPR vs. NJOY LEAPR\nUsing thin triangle approx.')
        ########################################################################
        # Compar S(a,b) errors
        ########################################################################
        plotErrorBetaForVariousAlpha(alphas,betas,MY_CONTIN2_SAB,NJOY_CONTIN2_SAB,A0,E,kbT,scalarMap)
        finishPlotting(colorBar,'S(a,b) error for H in H2O, generated using my LEAPR vs. NJOY LEAPR\nUsing thin triangle approx.')
        ########################################################################


    if option == 3:
        ########################################################################
        # COMPARE MY LEAPR DELTA FUNCS AGAINST MY LEAPR CONTIN2 FUNCS
        ########################################################################
        MY_DELTA_SAB   = getSAB(alphas,betas,continRho,NJOY_LEAPR=False,\
                                fullRedo=fullRedo,width=None,oscE=oscE,oscW=oscW)
        MY_CONTIN2_SAB = getSAB(alphas,betas,continRho,NJOY_LEAPR=False,\
                                fullRedo=fullRedo,width=2,oscE=oscE,oscW=oscW)
        ########################################################################
        # Compar S(a,b) values
        plotBetaForVariousAlpha(alphas,betas,MY_DELTA_SAB,A0,E,kbT,scalarMap,'.',False)
        plotBetaForVariousAlpha(alphas,betas,MY_CONTIN2_SAB,A0,E,kbT,scalarMap,'.',True)
        finishPlotting(colorBar,'S(a,b) for H in H2O, generated using my LEAPR v\nUsing thin triangle approx. vs delta func. approx.')
        ########################################################################
        # Compar S(a,b) errors
        ########################################################################
        plotErrorBetaForVariousAlpha(alphas,betas,MY_DELTA_SAB,MY_CONTIN2_SAB,A0,E,kbT,scalarMap)
        finishPlotting(colorBar,'S(a,b) error for H in H2O, generated using my LEAPR v\nUsing thin triangle approx. vs delta func. approx.')
        ########################################################################
    

    if option == 4:
        ########################################################################
        # COMPARE NJOY LEAPR DELTA FUNCS AGAINST NJOY LEAPR CONTIN2 FUNCS
        ########################################################################
        NJOY_DELTA_SAB   = getSAB(alphas,betas,continRho,NJOY_LEAPR=True,\
                                fullRedo=fullRedo,width=None,oscE=oscE,oscW=oscW)
        NJOY_CONTIN2_SAB = getSAB(alphas,betas,continRho,NJOY_LEAPR=True,\
                                fullRedo=fullRedo,width=2,oscE=oscE,oscW=oscW)
        ########################################################################
        # Compar S(a,b) values
        plotBetaForVariousAlpha(alphas,betas,NJOY_DELTA_SAB,A0,E,kbT,scalarMap,'.',False)
        plotBetaForVariousAlpha(alphas,betas,NJOY_CONTIN2_SAB,A0,E,kbT,scalarMap,'.',True)
        finishPlotting(colorBar,'S(a,b) for H in H2O, generated using NJOY LEAPR v\nUsing thin triangle approx. vs delta func. approx.')
        ########################################################################
        # Compar S(a,b) errors
        ########################################################################
        plotErrorBetaForVariousAlpha(alphas,betas,NJOY_DELTA_SAB,NJOY_CONTIN2_SAB,A0,E,kbT,scalarMap)
        finishPlotting(colorBar,'S(a,b) error for H in H2O, generated using NJOY LEAPR v\nUsing thin triangle approx. vs delta func. approx.')
        ###########################################################################








if __name__=="__main__":
    question = \
    " (1)  My  LEAPR delta function vs NJOY LEAPR delta function \n" + \
    " (2)  My  LEAPR thin triangle  vs NJOY LEAPR thin triangle \n"  + \
    " (3)  My  LEAPR delta function vs  MY  LEAPR thin triangle \n"  + \
    " (4) NJOY LEAPR delta function vs NJOY LEAPR thin triangle \n"  + \
    "What option do you want?: "


    option = int(input(question))
    runAlphasBetas(option)
