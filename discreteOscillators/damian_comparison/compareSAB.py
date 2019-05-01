import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx



def prepPlot(alpha):
    cnorm = colors.Normalize(vmin=0,vmax=1*len(alpha)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) 
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('gnuplot2')) 
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('Accent_r')) 
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('Set1')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alpha))])
    colorBar = plt.contourf([[0,0],[0,0]], alpha, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    plt.yscale('log'); plt.show()



with open("test9_sab") as f:
    test9_sab = [float(x) for x in f.readline().split()]
    test9_alpha = [float(x) for x in f.readline().split()]
    test9_beta = [float(x) for x in f.readline().split()]
with open("damianLEAPR_sab") as f:
    damian_sab = [float(x) for x in f.readline().split()]
    damian_alpha = [float(x) for x in f.readline().split()]
    damian_beta = [float(x) for x in f.readline().split()]

for i in range(len(test9_alpha)):
    assert(abs(test9_alpha[i]-damian_alpha[i])<1e-6)
for i in range(len(test9_beta)):
    assert(abs(test9_beta[i]-damian_beta[i])<1e-6)

alpha = test9_alpha
beta = test9_beta


scalarMap, colorBar = prepPlot(alpha)


fig = plt.figure()
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4], xticklabels=[],ylim=(2e-9,3e1))
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4],ylim=(2e-9,3e1))
ax1.set_yscale('log')
ax2.set_yscale('log')

#plt.subplot(1,2,1)

for a in range(len(alpha)):
    sab_to_plot = [test9_sab[b+a*len(beta)] for b in range(len(beta))]
    ax1.plot(beta,sab_to_plot,color=scalarMap.to_rgba(a),alpha=1.0)




#plt.subplot(1,2,2)
for a in range(len(alpha)):
    sab_to_plot = [damian_sab[b+a*len(beta)] for b in range(len(beta))]
    ax2.plot(beta,sab_to_plot,color=scalarMap.to_rgba(a),alpha=1.0)



finishPlotting(colorBar)
#plt.show()

#plt.legend(loc='best')
#plt.yscale('log')
#plt.show()





