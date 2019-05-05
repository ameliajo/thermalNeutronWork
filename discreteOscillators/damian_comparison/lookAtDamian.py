import matplotlib.pyplot as plt
import numpy as np


rawX = []
rawY = []

with open("ane_damian_data.txt") as f:
    vals = f.readlines()
    for val in vals:
        valVec = val.split()
        if len(valVec) == 0: break
        rawX.append(float(valVec[0]))
        rawY.append(float(valVec[1]))

rawX2 = []
rawY2 = []

with open("ane_damian2_data.txt") as f:
    vals = f.readlines()
    for val in vals:
        valVec = val.split()
        if len(valVec) == 0: break
        rawX2.append(float(valVec[0]))
        rawY2.append(float(valVec[1]))




def interpolate(X,Y,x):
    for i in range(len(X)-1):
        if X[i] <= x and x < X[i+1]:
            m = (Y[i+1]-Y[i])/(X[i+1]-X[i])
            b = Y[i]
            return m*(x-X[i])+b
    return 0.0

X = np.linspace(0.0,rawX[-1],200)
Y  = [0.0]*len(X)
Y2 = [0.0]*len(X)
for i in range(len(X)):
    Y[i]  = interpolate(rawX,rawY,X[i])
    Y2[i] = interpolate(rawX2,rawY2,X[i])

invArea = 1.0/np.trapz(rawY,x=rawX)
rawY = [x*invArea for x in rawY]
invArea = 1.0/np.trapz(Y,x=X)
Y = [x*invArea for x in Y]

invArea = 0.522421/np.trapz(rawY2,x=rawX2)
rawY2 = [x*invArea for x in rawY2]
invArea = 0.522421/np.trapz(Y2,x=X)
Y2 = [x*invArea for x in Y2]

#plt.plot(rawX,rawY)
#plt.plot(X,Y)
#plt.plot(rawX2,rawY2)
#plt.plot(X,Y2)
#plt.show()
#exit()
#inv_area = 1.0/np.trapz(Y,x=X)
#Y = [val*inv_area for val in Y]
#print(np.trapz(Y,x=X))



endf_X = np.linspace(0.00255,66*0.00255,67)
endf_Y = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02, .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091, .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542, .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804, .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535, .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, .04244, .042, 0.0]
inv_area = 0.444444/np.trapz(endf_Y,x=endf_X)
endf_Y = [val*inv_area for val in endf_Y]
#print(np.trapz(endf_Y,x=endf_X))


plt.rcParams.update({'font.size': 12})
plt.plot(X,Y,label='Full Model')
#plt.plot(endf_X,endf_Y,'r',label='LEAPR Test #9')
plt.plot(X,Y2,'r',label='Reduced Model')
#plt.plot([0.205,0.205],[0.0,6.0],'r')
#plt.plot([0.480,0.480],[0.0,6.0],'r')
plt.plot([0.205,0.205],[0.0,6.0],'r')
plt.plot([0.430,0.430],[0.0,6.0],'r')
plt.xlabel('Energy (eV)')
plt.ylabel('Phonon Distribution (arbitrary units)')
plt.legend(loc='upper right')
plt.show()


print(Y2)


