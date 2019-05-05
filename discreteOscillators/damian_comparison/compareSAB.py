import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np



def prepPlot(alpha):
    cnorm = colors.Normalize(vmin=0,vmax=alpha[-1]+5)# tab20 gnuplot2 Set1
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('Accent_r')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alpha))])
    colorBar = plt.contourf([[0,0],[0,0]], alpha, cmap=mymap,norm=cnorm)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar,ticks=np.arange(0,max(alpha),5),extend='both').ax.set_ylabel('alpha values')
    plt.yscale('log'); 
    plt.show()



#with open("test9_sab") as f:
with open("damianLEAPR2_sab") as f:
    test9_sab   = [float(x) for x in f.readline().split()]
    test9_alpha = [float(x) for x in f.readline().split()]
    test9_beta  = [float(x) for x in f.readline().split()]
with open("damianLEAPR_sab") as f:
    damian_sab   = [float(x) for x in f.readline().split()]
    damian_alpha = [float(x) for x in f.readline().split()]
    damian_beta  = [float(x) for x in f.readline().split()]

for i in range(len(test9_alpha)):
    assert(abs(test9_alpha[i]-damian_alpha[i])<1e-6)
for i in range(len(test9_beta)):
    assert(abs(test9_beta[i]-damian_beta[i])<1e-6)

alpha = test9_alpha
beta  = test9_beta

scalarMap, colorBar = prepPlot(np.linspace(alpha[0],alpha[-1],len(alpha)))

fig = plt.figure()
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4], xticklabels=[],ylim=(2e-9,1e3),ylabel='S(a,-b)',yscale='log')
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4],ylim=(2e-9,1e3),ylabel='S(a,-b)',yscale='log')


for a in range(len(alpha)):
    sab_to_plot = [test9_sab[b+a*len(beta)] for b in range(len(beta))]
    ax1.plot(beta,sab_to_plot,color=scalarMap.to_rgba(a),alpha=0.6)

for a in range(len(alpha)):
    sab_to_plot = [damian_sab[b+a*len(beta)] for b in range(len(beta))]
    ax2.plot(beta,sab_to_plot,color=scalarMap.to_rgba(a),alpha=0.6)

plt.xlabel('beta')
plt.show()


plt.rcParams.update({'font.size': 12})

for a in range(len(alpha)):
    #sab_to_plot = [test9_sab[b+a*len(beta)] for b in range(len(beta))]
    sab_to_plot = [damian_sab[b+a*len(beta)] for b in range(len(beta))]
    plt.plot(beta,sab_to_plot,color=scalarMap.to_rgba(a),alpha=0.6)
plt.xlabel('beta')
plt.ylabel('S(a,-b)')
plt.yscale('log')
plt.show()



