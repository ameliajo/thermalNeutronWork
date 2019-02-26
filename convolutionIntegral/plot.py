import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt



def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20c')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    #scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab10')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('iteration')
    #ax.set_facecolor('xkcd:light grey blue') # off white
    #ax.set_facecolor('xkcd:very light blue') # off white
    plt.yscale('log')
    plt.show()


      
