import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

def findArea(energy,rho):
    y_int = integrate.cumtrapz(rho, energy, initial=0)
    return y_int[-1]


def findAreaUnderDelta(energy,delta,delta_i_1,indexHalfWidth):
    return integrate.cumtrapz(                                       \
            delta[ delta_i_1-indexHalfWidth:delta_i_1+indexHalfWidth+1], \
            energy[delta_i_1-indexHalfWidth:delta_i_1+indexHalfWidth+1], initial=0)[-1]

def delta(index,delta_i_1,deltaWgt,spacing,order):
    h = (1/order)
    for distance in range(order):
        if abs(index-delta_i_1) == distance: 
            return ((order-distance)*h/order)*deltaWgt/spacing
    return 0.0


deltaFuncLocation1 = 0.205
deltaFuncLocation2 = 0.48
spacing = 0.005
delta_i_1 = int(deltaFuncLocation1/spacing)
delta_i_2 = int(deltaFuncLocation2/spacing)
endRho = 0.53
numSpaces = endRho/spacing+1
print(numSpaces)

# I want the area of my delta function to be 5
deltaWgt1 = 0.16667
deltaWgt2 = 0.33333
# I want the area of my continuous part to be 5
continWgt = 0.44444


waterRho = [ 0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 
  0.0165, 0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
  0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 0.1218, 
  0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 0.0839, 0.0807, 
  0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 0.06652, 0.065, 0.0634, 
  0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 0.05462, 0.0535, 0.0525, 0.0515, 
  0.05042, 0.04934, 0.04822, 0.04706, 0.0459, 0.04478, 0.04366, 0.04288, 
  0.04244, 0.042, 0.0 ]
waterEnergy = [0.00255 * i for i in range(len(waterRho))]


energy = np.linspace(0,endRho,numSpaces)
rho = [0.0]*int(numSpaces)
numTriangles = 6
#print(energy)
for triangleWidth in range(1,numTriangles):
    delta1Vec = [delta(i,delta_i_1,deltaWgt1,spacing,triangleWidth) for i in range(len(energy))]
    delta2Vec = [delta(i,delta_i_2,deltaWgt2,spacing,triangleWidth) for i in range(len(energy))]
    delta1FuncArea = findAreaUnderDelta(energy,delta1Vec,delta_i_1,triangleWidth)
    assert(abs(delta1FuncArea - deltaWgt1) < 1e-7)
    delta2FuncArea = findAreaUnderDelta(energy,delta2Vec,delta_i_2,triangleWidth)
    assert(abs(delta2FuncArea - deltaWgt2) < 1e-7)

    continVec = [waterRho[int(i*spacing/0.00255)] if i*spacing/0.00255 < len(waterRho) else 0.0 for i in range(len(energy))]
    continArea = findArea(energy,continVec)
    continVec = [entry*continWgt / (continArea) for entry in continVec]
    continArea = findArea(energy,continVec)
    assert(abs(continArea - continWgt) < 1e-7)
    
    for i in range(len(energy)): rho[i] = continVec[i]+delta1Vec[i]+delta2Vec[i]
    plt.plot(energy,rho,label="Base Width = "+str(2*triangleWidth*spacing))

    totalArea = findArea(energy,rho)
    
    #print(rho)






#findAreaUnderDistribution(energy,deltaVec)


plt.title("H in H2O distribution, with Delta Functions made Continuous")
plt.legend(loc="upper left")
plt.xlabel("Energy (eV)")
plt.ylabel("Freq. Distribution")
plt.show()
