import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np



def prepPlot(vec):
    plt.rcParams.update({'font.size': 12})
    cnorm = colors.Normalize(vmin=0,vmax=len(vec)+1)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('rainbow')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vec))])
    colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar

def finishPlotting(colorBar,log):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('alpha values')
    if log:
        plt.yscale('log')
    plt.show()

def getThisBadBoy(name,outputAll):
    with open('../outputs/'+name+'_output') as f:
        sab    = [float(x) for x in f.readline().split()]
        alphas = [float(x) for x in f.readline().split()]
        betas  = [float(x) for x in f.readline().split()]

    sab_ordered = [ [sab[a*len(betas)+b] for b in range(len(betas))] \
                                        for a in range(len(alphas))]
    if outputAll:
        return alphas,betas,sab_ordered 
    return sab_ordered


names = ["leapr_H_in_ZrH","H_in_ZrH_gamma","H_in_ZrH2_epsilon"]
alphas, betas, sab_leapr = getThisBadBoy(names[0],outputAll=True)
sab_ZrH  = getThisBadBoy(names[1],outputAll=False)
sab_ZrH2 = getThisBadBoy(names[2],outputAll=False)
nalpha = len(alphas)
rbeta = range(len(betas))

scalarMap,colorBar = prepPlot(alphas)

plot_neg_side_sab    = False
plot_neg_side_errors = False
plot_full_sab        = False
plot_xs              = True

if plot_neg_side_sab:
    for a in range(0,nalpha,10):
        plt.plot(betas,sab_leapr[a],color=scalarMap.to_rgba(a),linestyle="-")
        plt.plot(betas,sab_ZrH[a],color=scalarMap.to_rgba(a),linestyle=":")
        plt.plot(betas,sab_ZrH2[a],color=scalarMap.to_rgba(a),linestyle="-.")
    plt.xlabel('beta')
    plt.ylabel('S(a,-b)')
    finishPlotting(colorBar,log=True)

 

if plot_neg_side_errors:
    for a in range(0,nalpha,10):
        diff_ZrH2_ZrH = [sab_ZrH2[a][b]-sab_ZrH[a][b] for b in rbeta]
        #pct_diff_ZrH2_ZrH = [100*(sab_ZrH2[a][b]-sab_ZrH[a][b])/sab_ZrH2[a][b] for b in range(nbeta)]
        plt.plot(betas,diff_ZrH2_ZrH,color=scalarMap.to_rgba(a))
    plt.xlabel('beta')
    plt.ylabel('S(a,-b) error')
    finishPlotting(colorBar,log=False)


if plot_full_sab:
    full_betas = [-x for x in betas[::-1][:-1]] + betas
    for a in range(0,nalpha,10):
        endf = sab_leapr[a][::-1][:-1] + [np.exp(-betas[b])*sab_leapr[a][b] for b in rbeta]
        ZrH  = sab_ZrH[a][::-1][:-1]   + [np.exp(-betas[b])*sab_ZrH[a][b]   for b in rbeta]
        ZrH2 = sab_ZrH2[a][::-1][:-1]  + [np.exp(-betas[b])*sab_ZrH2[a][b]  for b in rbeta]

        plt.plot(full_betas,endf,color=scalarMap.to_rgba(a),linestyle="-")
        plt.plot(full_betas,ZrH,color=scalarMap.to_rgba(a),linestyle=":")
        plt.plot(full_betas,ZrH2,color=scalarMap.to_rgba(a),linestyle="-.")

    plt.xlabel('beta')
    plt.ylabel('S(a,-b)')
    finishPlotting(colorBar,log=True)


def getIndex(val,vec):
    for i in range(len(vec)-1):
        if vec[i] <= val <= vec[i+1]:
            return i, (val-vec[i])/(vec[i+1]-vec[i])
    return len(vec)-1,0


def getSABval(alpha,beta,alphas,betas,sab):
    a,aFrac = getIndex(alpha,alphas)
    b,bFrac = getIndex(abs(beta),betas)

    aL_bL, aR_bL = sab[a][b], sab[a][b]
    aL_bR, aR_bR = sab[a][b], sab[a][b]

    if a < len(alphas)-1:
        aR_bL = sab[a+1][b]
        aR_bR = sab[a+1][b]
    if b < len(betas)-1:
        aL_bR = sab[a][b+1]
        aR_bR = sab[a][b+1]
    if a < len(alphas)-1 and b < len(betas)-1:
        aR_bR = sab[a+1][b+1]

    one = aL_bL*(1-aFrac) + aR_bL*aFrac
    two = aL_bR*(1-aFrac) + aR_bR*aFrac

    three = one*(1-bFrac)+two*bFrac

    if beta < 0:
        return three
    return three*np.exp(-beta)







if plot_xs:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6
    E = 0.5
    mu = 0.5
    full_betas = [-x for x in betas[::-1][:-1]] + betas
    for a in range(0,nalpha,10):
        endf = [np.exp(-betas[b])*sab_leapr[a][b]  for b in rbeta][::-1][:-1] + sab_leapr[a]
        ZrH  = [np.exp(-betas[b])*sab_ZrH[a][b]  for b in rbeta][::-1][:-1] + sab_ZrH[a]
        ZrH2 = [np.exp(-betas[b])*sab_ZrH2[a][b] for b in rbeta][::-1][:-1] + sab_ZrH2[a]

    Ep_grid = np.linspace(0,1,300)
    leapr_sigma = []
    ZrH1_sigma = []
    ZrH2_sigma = []
    for Ep in Ep_grid:
        beta = (Ep-E)/kbT
        alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
        leapr_sigma.append(getSABval(alpha,beta,alphas,betas,sab_leapr))
        ZrH1_sigma.append(getSABval(alpha,beta,alphas,betas,sab_ZrH))
        ZrH2_sigma.append(getSABval(alpha,beta,alphas,betas,sab_ZrH2))

    plt.plot(Ep_grid,leapr_sigma,label='leapr')
    plt.plot(Ep_grid,ZrH1_sigma,label='ZrH')
    plt.plot(Ep_grid,ZrH2_sigma,label='ZrH2')
    plt.show()






    #plt.xlabel('beta')
    #plt.ylabel('S(a,-b)')
    #finishPlotting(colorBar,log=True)



    









