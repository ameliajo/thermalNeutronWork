import matplotlib.pyplot as plt
import numpy as np

def delta(n,deltaFuncLocation):
    print(abs(n - deltaFuncLocation) < 1e-7)
    return 1.0 if abs(n - deltaFuncLocation) < 1e-7 else 0.0

def delta2(n,deltaFuncLocation):
    return 1.0 if abs(n - deltaFuncLocation) < 1e-7 else \
           0.5 if abs(n - abs(deltaFuncLocation-spacing)) < 1e-7 else \
           0.5 if abs(n - abs(deltaFuncLocation+spacing)) < 1e-7 \
           else 0.0

def delta3(n,deltaFuncLocation):
    return 1.0 if abs(n - deltaFuncLocation) < 1e-7 else \
           2/3 if abs(n - abs(deltaFuncLocation-1*spacing)) < 1e-7 else \
           1/3 if abs(n - abs(deltaFuncLocation-2*spacing)) < 1e-7 else \
           2/3 if abs(n - abs(deltaFuncLocation+1*spacing)) < 1e-7 else \
           1/3 if abs(n - abs(deltaFuncLocation+2*spacing)) < 1e-7 \
           else 0.0




deltaFuncLocation = 0.5
spacing = 0.01
endRho = 1.0
numSpaces = endRho/spacing

energy = np.linspace(0,endRho-spacing,numSpaces)
deltaVec = [1.0+delta(i,deltaFuncLocation) for i in energy]
plt.plot(energy,deltaVec)

deltaVec = [1.0+delta2(i,deltaFuncLocation) for i in energy]
plt.plot(energy,deltaVec)

deltaVec = [1.0+delta3(i,deltaFuncLocation) for i in energy]
plt.plot(energy,deltaVec)




plt.show()

