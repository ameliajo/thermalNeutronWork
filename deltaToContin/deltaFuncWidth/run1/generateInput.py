import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

def findArea(energy,rho):
    y_int = integrate.cumtrapz(rho, energy, initial=0)
    return y_int[-1]


def findAreaUnderDelta(energy,delta,delta_i,indexHalfWidth):
    return integrate.cumtrapz(                                       \
            delta[ delta_i-indexHalfWidth:delta_i+indexHalfWidth+1], \
            energy[delta_i-indexHalfWidth:delta_i+indexHalfWidth+1], initial=0)[-1]

def delta(index,delta_i,deltaWgt,spacing,order):
    h = (1/order)
    for distance in range(order):
        if abs(index-delta_i) == distance: 
            return ((order-distance)*h/order)*deltaWgt/spacing
    return 0.0


deltaFuncLocation = 0.5
spacing = 0.01
delta_i = int(deltaFuncLocation/spacing)
endRho = 1.0+spacing
numSpaces = endRho/spacing

# I want the area of my delta function to be 5
deltaWgt = 2.0
# I want the area of my continuous part to be 5
continWgt = 1.0


energy = np.linspace(0,endRho-spacing,numSpaces)
rho = [0.0]*int(numSpaces)
numTriangles = 6
for triangleWidth in range(1,numTriangles):
    deltaVec = [delta(i,delta_i,deltaWgt,spacing,triangleWidth) for i in range(len(energy))]
    deltaFuncArea = findAreaUnderDelta(energy,deltaVec,delta_i,triangleWidth)
    assert(abs(deltaFuncArea - deltaWgt) < 1e-2)

    continVec = [continWgt/1.0 for i in range(len(energy))]
    continArea = findArea(energy,continVec)
    assert(abs(continArea - continWgt) < 1e-2)
    
    for i in range(len(energy)): rho[i] = continVec[i]+deltaVec[i]
    plt.plot(energy,rho,label="Base Width = "+str(2*triangleWidth*spacing))

    totalArea = findArea(energy,rho)
    
    #print(deltaFuncArea,continArea,totalArea)
    print(energy)
    print(rho)






#findAreaUnderDistribution(energy,deltaVec)



plt.title("Mock Freq. Dist. for Single Delta Func. with Small Contin")
plt.legend(loc="upper right")
plt.xlabel("Energy (eV)")
plt.ylabel("Freq. Distribution")
plt.show()
