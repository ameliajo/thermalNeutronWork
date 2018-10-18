import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

def findAreaUnderDistribution(energy,rho):
    y_int = integrate.cumtrapz(rho, energy, initial=0)
    print(y_int)
    plt.plot(energy,rho,'r')
    plt.plot(energy,y_int,'b')
    plt.show()


def findAreaUnderDelta(energy,rho,delta_i,indexHalfWidth):
    y_int = integrate.cumtrapz(rho, energy, initial=0)
    y_int = integrate.cumtrapz(rho[delta_i-indexHalfWidth:delta_i+indexHalfWidth+1], energy[delta_i-indexHalfWidth:delta_i+indexHalfWidth+1], initial=0)
    return y_int




def delta(index,delta_i,deltaWgt,spacing):
    return deltaWgt/spacing if index == delta_i else 0.0

def delta2(index,delta_i,deltaWgt,spacing):
    return 0.50*deltaWgt/spacing if index == delta_i else \
           0.25*deltaWgt/spacing if abs(index-delta_i) == 1 else 0.0

def delta3(index,delta_i,deltaWgt,spacing):
    h = (1/3)
    return h*deltaWgt/spacing if index == delta_i else \
           (2*h/3)*deltaWgt/spacing if abs(index-delta_i) == 1 else \
           (1*h/3)*deltaWgt/spacing if abs(index-delta_i) == 2 else 0.0

def delta4(index,delta_i,deltaWgt,spacing):
    h = (1/4)
    return h*deltaWgt/spacing if index == delta_i else \
           (3*h/4)*deltaWgt/spacing if abs(index-delta_i) == 1 else \
           (2*h/4)*deltaWgt/spacing if abs(index-delta_i) == 2 else \
           (1*h/4)*deltaWgt/spacing if abs(index-delta_i) == 3 else 0.0




deltaFuncLocation = 0.5
spacing = 0.05
delta_i = int(deltaFuncLocation/spacing)
endRho = 1.0+spacing
numSpaces = endRho/spacing

# I want the area of my delta function to be 5
deltaWgt = 2.0
# I want the area of my continuous part to be 5
continWgt = 1.0


energy = np.linspace(0,endRho-spacing,numSpaces)

deltaVec = [delta(i,delta_i,deltaWgt,spacing) for i in range(len(energy))]
plt.plot(energy,deltaVec)
deltaFuncArea = findAreaUnderDelta(energy,deltaVec,delta_i,1)
assert(abs(deltaFuncArea[-1] - deltaWgt) < 1e-15)

deltaVec = [delta2(i,delta_i,deltaWgt,spacing) for i in range(len(energy))]
plt.plot(energy,deltaVec)
deltaFuncArea = findAreaUnderDelta(energy,deltaVec,delta_i,2)
assert(abs(deltaFuncArea[-1] - deltaWgt) < 1e-15)

deltaVec = [delta3(i,delta_i,deltaWgt,spacing) for i in range(len(energy))]
plt.plot(energy,deltaVec)
deltaFuncArea = findAreaUnderDelta(energy,deltaVec,delta_i,3)
assert(abs(deltaFuncArea[-1] - deltaWgt) < 1e-15)

deltaVec = [delta4(i,delta_i,deltaWgt,spacing) for i in range(len(energy))]
plt.plot(energy,deltaVec)
deltaFuncArea = findAreaUnderDelta(energy,deltaVec,delta_i,4)
assert(abs(deltaFuncArea[-1] - deltaWgt) < 1e-15)




"""

findAreaUnderDistribution(energy,deltaVec)



"""
plt.show()
