import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from mpl_toolkits import mplot3d




def prepPlot(vec):
    plt.rcParams.update({'font.size': 12})
    cnorm = colors.Normalize(vmin=0,vmax=len(vec))
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('rainbow')) 
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(vec))])
    colorBar = plt.contourf([[0,0],[0,0]], vec, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar

def finishPlotting(colorBar,log,title):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel(title)
    if log: plt.yscale('log')
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

def getIndex(val,vec):
    for i in range(len(vec)-1):
        if vec[i] <= val <= vec[i+1]:
            return i, (val-vec[i])/(vec[i+1]-vec[i])
    return len(vec)-1,0


def getSABval(alpha,beta,alphas,betas,sab):
    a,aFrac = getIndex(alpha,alphas)
    b,bFrac = getIndex(abs(beta),betas)

    aL_bL, aR_bL, aL_bR, aR_bR = sab[a][b], sab[a][b], sab[a][b], sab[a][b]
    #aL_bL, aR_bL, aL_bR, aR_bR = 0.0,0.0,0.0,0.0

    if a < len(alphas)-1: aR_bL = sab[a+1][b]; aR_bR = sab[a+1][b]
    if b < len(betas)-1:  aL_bR = sab[a][b+1]; aR_bR = sab[a][b+1]
    if a < len(alphas)-1 and b < len(betas)-1: aR_bR = sab[a+1][b+1]

    sabVal = ( aL_bL*(1-aFrac) + aR_bL*aFrac )*(1-bFrac) + \
             ( aL_bR*(1-aFrac) + aR_bR*aFrac )*bFrac

    return sabVal if beta < 0 else sabVal*np.exp(-beta)






names = ["leapr_H_in_ZrH","H_in_ZrH_gamma","H_in_ZrH2_epsilon","Zr_in_ZrH2_epsilon"]
alphas, betas, sab_leapr = getThisBadBoy(names[0],outputAll=True)
sab_H_ZrH  = getThisBadBoy(names[1],outputAll=False)
sab_H_ZrH2 = getThisBadBoy(names[2],outputAll=False)
Zr_alphas, Zr_betas, sab_Zr_ZrH2 = getThisBadBoy(names[3],outputAll=True)
nalpha = len(alphas)
rbeta = range(len(betas))

scalarMap,colorBar = prepPlot(alphas)

plot_neg_side_sab     = 0
plot_neg_side_errors  = 0
plot_full_sab         = 0
plot_xs_various_dists = 0
plot_xs_various_mu    = 0
plot_xs_various_dists_and_mu = 0
compareWithCouch = 1
compareWithPurohit = 0
plot_integrate_over_mu = 0




if plot_neg_side_sab:
    for a in range(0,nalpha,10):
        plt.plot(betas,sab_leapr[a],color=scalarMap.to_rgba(a),linestyle="-")
        plt.plot(betas,sab_H_ZrH[a],color=scalarMap.to_rgba(a),linestyle=":")
        plt.plot(betas,sab_H_ZrH2[a],color=scalarMap.to_rgba(a),linestyle="-.")
    plt.xlabel('beta'); plt.ylabel('S(a,-b)')
    finishPlotting(colorBar,log=True,title='alpha values')

 

if plot_neg_side_errors:
    for a in range(0,nalpha,10):
        diff_ZrH2_ZrH = [sab_H_ZrH2[a][b]-sab_H_ZrH[a][b] for b in rbeta]
        plt.plot(betas,diff_ZrH2_ZrH,color=scalarMap.to_rgba(a))
    plt.xlabel('beta'); plt.ylabel('S(a,-b) error')
    finishPlotting(colorBar,log=True,title='alpha values')


if plot_full_sab:
    full_betas = [-x for x in betas[::-1][:-1]] + betas
    for a in range(0,nalpha,10):
        endf = sab_leapr[a][::-1][:-1] + [np.exp(-betas[b])*sab_leapr[a][b] for b in rbeta]
        ZrH  = sab_H_ZrH[a][::-1][:-1]   + [np.exp(-betas[b])*sab_H_ZrH[a][b]   for b in rbeta]
        ZrH2 = sab_H_ZrH2[a][::-1][:-1]  + [np.exp(-betas[b])*sab_H_ZrH2[a][b]  for b in rbeta]
        plt.plot(full_betas,endf,color=scalarMap.to_rgba(a),linestyle="-")
        plt.plot(full_betas,ZrH,color=scalarMap.to_rgba(a),linestyle=":")
        plt.plot(full_betas,ZrH2,color=scalarMap.to_rgba(a),linestyle="-.")
    plt.xlabel('beta'); plt.ylabel('S(a,-b)')
    finishPlotting(colorBar,log=True,title='alpha values')


if plot_xs_various_dists:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6
    E = 0.5; mu = 0.6

    full_betas = [-x for x in betas[::-1][:-1]] + betas
    Ep_grid = np.linspace(0,1.4*E,3000)
    leapr_sigma, ZrH1_sigma, ZrH2_sigma = [], [], []
    for Ep in Ep_grid:
        beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
        leapr_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_leapr))
        #ZrH1_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_H_ZrH))
        #ZrH2_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_H_ZrH2))

    plt.rcParams.update({'font.size': 14})
    plt.plot(Ep_grid,leapr_sigma,label='leapr')
    #plt.plot(Ep_grid,ZrH1_sigma,label='ZrH')
    plt.plot([E+.137,E+.137],[0,max(leapr_sigma)],'g',alpha=0.5)
    plt.plot([E,E],[0,max(leapr_sigma)],'r',alpha=0.5)
    plt.plot([E-.137,E-.137],[0,max(leapr_sigma)],'g',alpha=0.5)
    plt.plot([E-2*.137,E-2*.137],[0,max(leapr_sigma)],'m',alpha=0.5)
    plt.plot([E-3*.137,E-3*.137],[0,max(leapr_sigma)],'y',alpha=1.0)
    plt.xlabel("Outgoing energy E' (eV)")
    plt.ylabel("xs (b)")
    #plt.legend(loc='best')
    plt.title("Inelastic xs of H in ZrH")
    plt.show()



if plot_xs_various_mu:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6
    E = 0.5
    mus = np.linspace(-1,1,20)
    scalarMap,colorBar = prepPlot(mus)
    full_betas = [-x for x in betas[::-1][:-1]] + betas
    for i,mu in enumerate(mus):
        Ep_grid = np.linspace(0,1,300)
        leapr_sigma = []
        ZrH1_sigma = []
        ZrH2_sigma = []
        for Ep in Ep_grid:
            beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
            leapr_sigma.append(getSABval(alpha,beta,alphas,betas,sab_leapr))
            ZrH1_sigma.append(getSABval(alpha,beta,alphas,betas,sab_H_ZrH))
            ZrH2_sigma.append(getSABval(alpha,beta,alphas,betas,sab_H_ZrH2))

        plt.plot(Ep_grid,leapr_sigma,label=str(mu),color=scalarMap.to_rgba(i))
        #plt.plot(Ep_grid,ZrH1_sigma,label=str(mu),color=scalarMap.to_rgba(i))
        #plt.plot(Ep_grid,ZrH2_sigma,label=str(mu),color=scalarMap.to_rgba(i))
    plt.xlabel('beta'); plt.ylabel('xs')
    finishPlotting(colorBar,log=False,title='mu')

    
if plot_xs_various_dists_and_mu:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6
    E = 0.5; mu = 0.5
    full_betas = [-x for x in betas[::-1][:-1]] + betas
    Ep_grid = np.linspace(0,1,3000)
    mus = np.linspace(-1,1,10)
    scalarMap,colorBar = prepPlot(mus)
    for i,mu in enumerate(mus):
        leapr_sigma, ZrH1_sigma, ZrH2_sigma = [], [], []
        for Ep in Ep_grid:
            beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
            leapr_sigma.append(getSABval(alpha,beta,alphas,betas,sab_leapr))
            ZrH1_sigma.append(getSABval(alpha,beta,alphas,betas,sab_H_ZrH))
            ZrH2_sigma.append(getSABval(alpha,beta,alphas,betas,sab_H_ZrH2))

        plt.plot(Ep_grid,leapr_sigma,color=scalarMap.to_rgba(i),linestyle='-')
        plt.plot(Ep_grid,ZrH1_sigma,color=scalarMap.to_rgba(i),linestyle='-.')
        plt.plot(Ep_grid,ZrH2_sigma,color=scalarMap.to_rgba(i),linestyle=':')
    finishPlotting(colorBar,log=False,title='mu')


if compareWithCouch:
    couch_x = [1.048e2, 1.062e2, 1.077e2, 1.091e2, 1.106e2, 1.121e2, 1.135e2, 
      1.151e2, 1.165e2, 1.179e2, 1.194e2, 1.208e2, 1.223e2, 1.238e2, 1.242e2, 
      1.252e2, 1.255e2, 1.265e2, 1.273e2, 1.279e2, 1.276e2, 1.283e2, 1.285e2, 
      1.288e2, 1.291e2, 1.291e2, 1.297e2, 1.299e2, 1.301e2, 1.305e2, 1.305e2, 
      1.308e2, 1.310e2, 1.314e2, 1.316e2, 1.318e2, 1.321e2, 1.324e2, 1.326e2, 
      1.329e2, 1.332e2, 1.336e2, 1.341e2, 1.344e2, 1.339e2, 1.346e2, 1.350e2, 
      1.355e2, 1.353e2, 1.356e2, 1.363e2, 1.368e2, 1.380e2, 1.376e2, 1.385e2, 
      1.392e2, 1.394e2, 1.401e2, 1.404e2, 1.409e2, 1.421e2, 1.427e2, 1.430e2, 
      1.434e2, 1.447e2, 1.455e2, 1.458e2, 1.461e2, 1.463e2, 1.468e2, 1.470e2, 
      1.475e2, 1.476e2, 1.478e2, 1.479e2, 1.480e2, 1.483e2, 1.484e2, 1.490e2, 
      1.495e2, 1.502e2, 1.508e2, 1.512e2, 1.521e2, 1.531e2, 1.541e2, 1.550e2, 
      1.560e2, 1.563e2, 1.569e2, 1.574e2, 1.579e2, 1.585e2, 1.591e2, 1.596e2, 
      1.600e2, 1.611e2, 1.607e2, 1.623e2, 1.634e2, 1.630e2, 1.643e2, 1.650e2]

    couch_y = [4.193, 4.174, 4.112, 4.441, 5.067, 4.869, 4.751, 5.242, 5.832, 
      6.830, 7.759, 8.563, 9.706, 1.076e1, 1.167e1, 1.283e1, 1.4e1, 1.585e1, 
      1.74e1, 1.983e1, 1.827e1, 2.153e1, 2.273e1, 2.415e1, 2.651e1, 2.495e1, 
      2.838e1, 2.964e1, 3.089e1, 3.24e1, 3.365e1, 3.509e1, 3.638e1, 3.773e1, 
      3.894e1, 4.035e1, 4.164e1, 4.297e1, 4.428e1, 4.562e1, 4.71e1, 4.878e1, 
      5.24e1, 5.361e1, 4.971e1, 5.474e1, 5.594e1, 5.758e1, 5.671e1, 5.952e1,
      6.112e1, 6.185e1, 6.021e1, 6.124e1, 5.83e1, 5.685e1, 5.527e1, 5.357e1, 
      5.184e1, 5.078e1, 5.065e1, 5.2e1, 5.307e1, 5.427e1, 5.423e1, 5.283e1, 
      5.128e1, 5.015e1, 4.896e1, 4.73e1, 4.553e1, 4.412e1, 4.299e1, 4.182e1, 
      4.048e1, 3.888e1, 3.754e1, 3.647e1, 3.455e1, 3.288e1, 3.159e1, 3.026e1, 
      2.929e1, 2.766e1, 2.656e1, 2.559e1, 2.442e1, 2.269e1, 2.148e1, 2.035e1, 
      1.916e1, 1.79e1, 1.669e1, 1.513e1, 1.378e1, 1.219e1, 9.441, 1.103e1, 
      7.297, 4.042, 5.690, 3.039, 1.778]

    couch_x = [104.7594, 106.2107, 107.673, 109.137, 110.6031, 112.0647, 113.5268, 115.0588, 116.4963, 117.9495, 119.3660, 120.8451, 122.2789, 123.7947, 124.1664, 125.1530, 125.5276, 126.4999, 127.2590, 127.580, 127.8810, 128.3258, 128.519, 128.7605, 129.0878, 129.122, 129.6808, 129.8769, 130.129, 130.4696, 130.4885, 130.7533, 130.998, 131.418, 131.6211, 131.8237, 132.0831, 132.3682, 132.639, 132.9125, 133.230, 133.5724, 133.9451, 134.108, 134.3820, 134.631, 134.9893, 135.2562, 135.4711, 135.6278, 136.2645, 136.841, 137.553, 137.997, 138.4538, 139.2236, 139.3666, 140.1105, 140.4276, 140.933, 142.055, 142.662, 143.016, 143.4398, 144.7226, 145.509, 145.814, 146.0663, 146.339, 146.7605, 147.0067, 147.5072, 147.5790, 147.7668, 147.892, 148.022, 148.2977, 148.4372, 149.0468, 149.50, 150.243, 150.8262, 151.241, 152.1225, 153.0552, 154.1388, 155.029, 155.9964, 156.2691, 156.885, 157.431, 157.919, 158.5222, 159.0937, 159.5666, 160.004, 160.7104, 161.075, 162.2593, 163.0176, 163.4262, 164.2994, 165.0322]

    couch_y = [4.19250, 4.17375, 4.11204, 4.44097, 5.06729, 4.86935, 4.75083, 5.24170, 5.83198, 6.83033, 7.75886, 8.56287, 9.70643, 10.75927, 11.67330, 12.8271, 14.0033, 15.85457, 17.39978, 18.27, 19.82734, 21.5287, 22.72667, 24.1458, 24.94605, 26.51071, 28.37679, 29.6421, 30.88758, 32.4028, 33.6492, 35.0862, 36.3774, 37.7300, 38.9406, 40.3477, 41.6431, 42.9659, 44.2796, 45.6207, 47.0986, 48.7779, 49.7072, 52.40267, 53.6051, 54.7388, 55.9374, 56.7142, 57.5838, 59.5213, 61.1182, 61.85335, 61.2432, 60.2093, 58.3017, 56.84800, 55.27048, 53.5676, 51.8366, 50.7796, 50.65231, 51.999, 53.0661, 54.2746, 54.2294, 52.82550, 51.28184, 50.15232, 48.95780, 47.3022, 45.53440, 44.12203, 42.986, 41.8171, 40.482, 38.8794, 37.54442, 36.4735, 34.54669, 32.87747, 31.5876, 30.257, 29.29048, 27.65666, 26.56006, 25.58945, 24.4240, 22.6946, 21.48331, 20.3465, 19.1617, 17.89912, 16.69144, 15.12689, 13.77785, 12.18848, 11.02844, 9.4412, 7.29729, 5.690073, 4.04177, 3.039135, 1.778274]

    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6

    E = 0.1715
    mu =-0.515038 # 121 
    mu = 0.96814764037 # 14.5 

    full_betas = [-x for x in betas[::-1][:-1]] + betas
    Ep_grid = np.linspace(0,1.5*E,500)
    H_ZrH2_sigma, Zr_ZrH2_sigma  = [], []
    leapr_sigma = []
    for Ep in Ep_grid:
        beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
        H_ZrH2_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_H_ZrH2))
        leapr_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_leapr))

        #beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(90.436*kbT)
        #Zr_ZrH2_sigma.append(6.20/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,Zr_alphas,Zr_betas,sab_Zr_ZrH2))

    plt.rcParams.update({'font.size': 14})
    plt.plot([1000*(E-Ep) for Ep in Ep_grid],H_ZrH2_sigma,label='DFT')
    plt.plot([1000*(E-Ep) for Ep in Ep_grid],leapr_sigma,label='LEAPR')
    #plt.plot([1000*(E-Ep) for Ep in Ep_grid],Zr_ZrH2_sigma,label='Zr in ZrH2')
    plt.title('Comparing H in ZrH2 scattering distributions')

    plt.plot(couch_x,couch_y,label='Couch')

    plt.xlim([110,1000*E])
    plt.xlabel("Energy lost by neutron (meV)")
    plt.ylabel("xs (b)")
    plt.legend(loc='best')
    plt.show()




if compareWithPurohit:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6

    E = 0.527
    mu = 0.90630778703 # 25 
    mu = 0.76604444311 # 40

    full_betas = [-x for x in betas[::-1][:-1]] + betas
    Ep_grid = np.linspace(0.9*E,1.1*E,5)
    Ep_grid = np.linspace(0,1.5*E,500)
    H_ZrH2_sigma, Zr_ZrH2_sigma, totalSigma  = [], [], []
    for Ep in Ep_grid:
        beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
        H_ZrH2_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,alphas,betas,sab_H_ZrH2))
        beta = (Ep-E)/kbT; alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(90.436*kbT)
        Zr_ZrH2_sigma.append(6.20/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,Zr_alphas,Zr_betas,sab_Zr_ZrH2))
        #print(Ep,getSABval(alpha,beta,Zr_alphas,Zr_betas,sab_Zr_ZrH2))

        totalSigma.append(H_ZrH2_sigma[-1]+Zr_ZrH2_sigma[-1])
        #Zr_ZrH2_sigma.append(sigma_b/(2.0*kbT)*(Ep/E)**0.5*getSABval(alpha,beta,Zr_alphas,Zr_betas,sab_Zr_ZrH2))

    plt.plot([Ep for Ep in Ep_grid],H_ZrH2_sigma,label='H in ZrH2')
    #plt.plot([Ep for Ep in Ep_grid],Zr_ZrH2_sigma,label='Zr in ZrH2')
    #plt.plot([Ep for Ep in Ep_grid],totalSigma,label='Total ZrH2')
    #plt.plot(Ep_grid,Zr_ZrH2_sigma,label='Zr in ZrH2')

    #plt.plot([E,E],[0,max(totalSigma)],'r')

    #plt.xlim([0,100*E])
    plt.xlabel("Outgoing energy E' (meV)")
    plt.ylabel("xs (b)")
    plt.legend(loc='best')
    plt.show()






if plot_integrate_over_mu:
    sigma_b = 20.478
    A = 0.99917
    kbT = 8.617e-5*293.6
    E = 1.5; 
    muVals = np.linspace(0,1.0,5)

    full_betas = [-x for x in betas[::-1][:-1]] + betas
    E_grid = np.linspace(0.1,5,5)
    for E in E_grid:
        Ep_grid = np.linspace(0,8,500)
        ZrH1_sigma = []
        for Ep in Ep_grid:
            ZrH1contrib = []
            beta = (Ep-E)/kbT
            for j,mu in enumerate(muVals):
                alpha = (Ep+E-2*mu*(Ep*E)**0.5)/(A*kbT)
                xsTerms = sigma_b/(2.0*kbT)*(Ep/E)**0.5
                ZrH1contrib.append(xsTerms*getSABval(alpha,beta,alphas,betas,sab_H_ZrH))
            ZrH1_sigma.append(np.trapz(ZrH1contrib,x=muVals))
        #plt.plot([Ep/E for Ep in Ep_grid],ZrH1_sigma,label=str(E))
        plt.plot(Ep_grid,ZrH1_sigma,label=str(E))
    plt.xlabel("Outgoing energy E' (eV)")
    #plt.xlim([0,1.5])
    plt.ylabel("xs (b)")
    plt.legend(loc='best')
    plt.show()




























