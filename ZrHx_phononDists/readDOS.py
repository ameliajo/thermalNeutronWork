import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np



endfDist_H = [0., .000875, .0035, .008, .015, .0235, .0340, .046, .061, .078, .094, .116, .144, .1606, .1969, .2606, .3479, .3559, .3500, .3322, .3328, .2911, .1617, .1431, .1248, .09738, .06067, .1221, .1495, .07219, .01443, .0001, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .0499, 2.010, 3.560, 4.790, 5.995, 7.250, 8.550, 9.640, 11.91, 13.52, 16.04, 19.79, 26.10, 29.39, 30.82, 32.21, 31.75, 33.14, 35.65, 33.34, 36.27, 38.18, 38.75, 39.48, 28.99, 23.29, 25.18, 26.59, 27.86, 27.89, 29.44, 25.86, 23.33, 24.66, 27.51, 37.94, 60.77, 26.66, 18.54, 14.51, 11.48, 9.53, 7.53, 5.449, 3.838, 8.497, 0.]
endfDist_Zr = [0.,.08,.32,.70,1.40,2.15,3.10,4.50,6.30,8.40,11.0,14.2,16.4,21.3, \
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
endfGrid = [i for i in range(len(endfDist_H))]
invArea = 1.0/np.trapz(endfDist_H,x=endfGrid)
endfDist_H = [invArea*x for x in endfDist_H]
invArea = 1.0/np.trapz(endfDist_Zr,x=endfGrid)
endfDist_Zr = [invArea*x for x in endfDist_Zr]





def interpolate(X,Y,x):
    for i in range(len(X)-1):
        if X[i] <= x <= X[i+1]:
            m = (Y[i+1]-Y[i])/(X[i+1]-X[i])
            b = Y[i]
            return m*(x-X[i])+b
    return None

def writeDataToFile(outputFileName,uniform_X,uniform_Y):
    with open(outputFileName,'w') as f:
        f.write('Spacing = '+str(1e-3*(uniform_X[1]-uniform_X[0]))+'\n\n')
        f.write('Length  = '+str(len(uniform_X))+'\n\n')
        for i in range(0,len(uniform_Y)-6,6):
            for j in range(6): f.write(str('%.4E'%uniform_Y[i+j])+" ")
            f.write('\n')



filenames = ["Zr_alpha", "ZrH2_epsilon", "ZrH_gamma"]
filenames = ["ZrH2_epsilon", "ZrH_gamma"]

num = 200

cnorm = colors.Normalize(vmin=0,vmax=len(filenames)+1)
scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('rainbow')) 

for index,filename in enumerate(filenames):
    print("Processing",filename)
    x, y = [], []
    counter = 0; split = None
    with open('rawDOS/'+filename+'.txt') as f:
        lines = f.readlines()
        for line in lines:
            data = line.split()
            if len(x) > 1 and float(data[1])-x[-1] > 50: split = counter 
            counter += 1
            y.append(float(data[0]))
            x.append(float(data[1]))

    # Do the hydrogen piece
    uniform_X = np.linspace(x[0],x[split-1],num)
    uniform_Y = [interpolate(x,y,uniform_X[i]) if i > split else 0.0 \
                 for i in range(len(uniform_X))]
    inv_area = 1.0/np.trapz(uniform_Y,x=uniform_X)
    uniform_Y = [val*inv_area for val in uniform_Y]
    writeDataToFile('uniformDOS/'+filename+'_uniform_DOS_Zr.txt',uniform_X,uniform_Y)

    plt.plot(uniform_X,uniform_Y,color=scalarMap.to_rgba(index))

    # Do the zirconium piece
    uniform_X = np.linspace(x[0],x[-1],num)
    uniform_Y = [interpolate(x,y,uniform_X[i]) if i > split else 0.0 \
                 for i in range(len(uniform_X))]
    inv_area = 1.0/np.trapz(uniform_Y,x=uniform_X)
    uniform_Y = [val*inv_area for val in uniform_Y]
    writeDataToFile('uniformDOS/'+filename+'_uniform_DOS_H.txt',uniform_X,uniform_Y)
    
    plt.plot(uniform_X,uniform_Y,label=filename,color=scalarMap.to_rgba(index))

plt.plot(endfGrid,endfDist_H,label='ENDF',color=scalarMap.to_rgba(index+1))
plt.plot(endfGrid,endfDist_Zr,color=scalarMap.to_rgba(index+1))
plt.legend(loc='best')
plt.xlabel('Energy (meV)')
plt.ylabel('Phonon DOS (arbitrary)')
plt.show()
















