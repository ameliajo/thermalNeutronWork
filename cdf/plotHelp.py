import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+3)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('hot')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar,title,xlabel,ylabel):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('beta values')
    plt.title(title)
    ax.set_facecolor('xkcd:light grey blue') # off white
    #plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


