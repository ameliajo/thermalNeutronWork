import matplotlib.pyplot as plt
import numpy as np


endfDist = [0., .000875, .0035, .008, .015, .0235, .0340, .046, .061, .078, .094, .116, .144, .1606, .1969, .2606, .3479, .3559, .3500, .3322, .3328, .2911, .1617, .1431, .1248, .09738, .06067, .1221, .1495, .07219, .01443, .0001, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .0499, 2.010, 3.560, 4.790, 5.995, 7.250, 8.550, 9.640, 11.91, 13.52, 16.04, 19.79, 26.10, 29.39, 30.82, 32.21, 31.75, 33.14, 35.65, 33.34, 36.27, 38.18, 38.75, 39.48, 28.99, 23.29, 25.18, 26.59, 27.86, 27.89, 29.44, 25.86, 23.33, 24.66, 27.51, 37.94, 60.77, 26.66, 18.54, 14.51, 11.48, 9.53, 7.53, 5.449, 3.838, 8.497, 0.]
endfGrid = [i for i in range(len(endfDist))]
#invArea = 1.0/np.trapz(endfDist,x=endfGrid)
endfDist = [0.02*x for x in endfDist]

zr_in_zrh = [0.,.08,.32,.70,1.40,2.15,3.10,4.50,6.30,8.40,11.0,14.2,16.4,21.3, \
             27.97,39.79,57.36,63.10,67.14,69.42,76.32,73.70,43.53,42.53,41.67,\
             37.72,29.24,66.07,94.47,58.62,19.57,1.031,0.,0.,0.,0.,0.,0.,0.,0.,\
             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\
             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\
             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,\
             0.,0.,0.,0.,0.,0.,0.,0.,.00218,.0828,.1340,.1650,.1860,.2050,     \
             .2150,.2170,.2360,.2340,.239,.242,.235,.223,.221,.214,.198,.189,  \
             .170,.206,.207,.202,.193,.197,.200,.212,.225,.243,.253,.253,.262, \
             .228,.195,.194,.204,.263,.390,.170,.108,.0775,.0592,.0406,.0235,  \
             .0112,.00424,.000266,0.]
endfDist = [0.02*x for x in endfDist]
#zr_in_zrh_energy = [i*0.001 for i in range(161)]






def interpolate(X,Y,x):
    for i in range(len(X)-1):
        if X[i] <= x <= X[i+1]:
            m = (Y[i+1]-Y[i])/(X[i+1]-X[i])
            b = Y[i]
            return m*(x-X[i])+b
    return None

filenames = ["Zr_alpha", "ZrH_gamma", "ZrH2_epsilon"]
filenames = ["ZrH_gamma", "ZrH2_epsilon"]
names = ["ZrH    (DFT)","ZrH2  (DFT)","ZrH  (ENDF)"]


colors = ["skyblue",'#ff9292','#f2ee97','#afe5ad','#83b8f4','#ecb3d2']
colors = ['#83b8f4','#ff9292','#ffb366','#ecb3d2','#f2ee97','#afe5ad']

plt.rcParams.update({'font.size': 14})

for i,filename in enumerate(filenames):
    x, y = [], []
    counter = 0
    split = None
    with open(filename+'.txt') as f:
        lines = f.readlines()
        for line in lines:
            data = line.split()
            if len(x) > 1 and float(data[1])-x[-1] > 50:
                split = counter 

            counter += 1
            y.append(float(data[0]))
            x.append(float(data[1]))

    invArea = 1.0/np.trapz(y,x=x)
    y = [val*invArea for val in y]
    plt.plot(x,y,label=names[i],color=colors[i])
    plt.fill_between( x, y, color=colors[i], alpha=0.6)


i+=1
invArea = 1.0/np.trapz(endfDist,x=endfGrid)
endfDist = [val*invArea for val in endfDist]
invArea = 1.0/np.trapz(zr_in_zrh ,x=endfGrid)
zr_in_zrh = [val*invArea for val in zr_in_zrh]
plt.plot(endfGrid[:int(len(endfGrid)/2)],zr_in_zrh[:int(len(endfGrid)/2)],colors[i])
plt.fill_between( endfGrid[:int(len(endfGrid)/2)], zr_in_zrh[:int(len(endfGrid)/2)], color=colors[i], alpha=0.3)


plt.plot(endfGrid[int(len(endfGrid)/2):],endfDist[int(len(endfGrid)/2):],colors[i],label=names[-1])
plt.fill_between( endfGrid[int(len(endfGrid)/2):], endfDist[int(len(endfGrid)/2):], color=colors[i], alpha=0.3)
#plt.plot(endfGrid,endfDist,label='ENDF')
plt.legend(loc='best')
plt.xlabel('Energy (meV)')
plt.ylabel('Phonon DOS (arb. units)')

plt.show()
