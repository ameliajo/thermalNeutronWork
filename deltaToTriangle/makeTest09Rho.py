import matplotlib.pyplot as plt 
from scipy import integrate
import numpy as np

def findArea(energy,rho):
    y_int = integrate.cumtrapz(rho, energy, initial=0)
    return y_int[-1]


spacing = 0.00255

E = [0.0]
for i in range(200):
    E.append(E[-1]+spacing)



continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
  .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
  .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
  .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
  .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
  .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
  .04244, .042, 0.0]


def getPhononDist(width,continRho):

    # Add in the continuous part that we copy from Test09 phonon distribution
    rho = [continRho[i] if i < len(continRho) else 0.0 for i in range(201)]

    # Normalize the continuous part
    contin_Weight  = 0.5
    areaUnderContin = findArea(E,rho)
    rho = [x*contin_Weight/areaUnderContin for x in rho]


    # Add in the delta function
    delta_Energies = [0.205,0.48]
    delta_Weights  = [0.166667,0.333333]

    # Figure out which indices to put my delta functions at
    i0 = int(delta_Energies[0]/spacing)
    i1 = int(delta_Energies[1]/spacing)
    # Check to make sure that they're reasonable
    #print(E[i0],E[i1])

    i = int(width/2)
    for j in range(i):
        rho[i0-j] = (i-j)*delta_Weights[0]/(i**2*spacing)
        rho[i0+j] = (i-j)*delta_Weights[0]/(i**2*spacing)
        rho[i1-j] = (i-j)*delta_Weights[1]/(i**2*spacing)
        rho[i1+j] = (i-j)*delta_Weights[1]/(i**2*spacing)

    return rho






if __name__=="__main__":
    for i in range(2,12,2):
        plt.plot(E,getPhononDist(i,continRho))
        print(i)
        #print(['%.5E'%x for x in getPhononDist(i,continRho)] )
        print(np.array(getPhononDist(i,continRho)))
        print()
    plt.show()




