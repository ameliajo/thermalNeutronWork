import matplotlib.pyplot as plt
import numpy as np


def interpolate(X,Y,x):
    for i in range(len(X)-1):
        if X[i] <= x <= X[i+1]:
            m = (Y[i+1]-Y[i])/(X[i+1]-X[i])
            b = Y[i]
            return m*(x-X[i])+b
    return None

filenames = ["Zr_alpha", "ZrH_gamma", "ZrH2_epsilon"]




for filename in filenames:
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

    plt.plot(x,y,label=filename)
plt.legend(loc='best')
plt.show()
"""
num = 500
uniform_X = np.linspace(x[0],x[-1],num)
uniform_Y = [0.0]*num

for i in range(len(uniform_X)):
    uniform_Y[i] = interpolate(x,y,uniform_X[i])

inv_area = 1.0/np.trapz(uniform_Y,x=uniform_X)
uniform_Y = [val*inv_area for val in uniform_Y]


with open('uniformDOS/'+filename+'_uniform_DOS.txt','w') as f:
    f.write('Spacing = '+str(uniform_X[1]-uniform_X[0])+'\n\n')
    f.write('Length  = '+str(len(uniform_X))+'\n\n')
    for i in range(0,len(uniform_Y)-6,6):
        for j in range(6):
            f.write(str('%.4E'%uniform_Y[i+j])+" ")
        f.write('\n')



num = 100
uniform_X = np.linspace(x[0],x[split-1],num)
uniform_Y = [0.0]*num

for i in range(len(uniform_X)):
    uniform_Y[i] = interpolate(x,y,uniform_X[i])

inv_area = 1.0/np.trapz(uniform_Y,x=uniform_X)
uniform_Y = [val*inv_area for val in uniform_Y]

with open('uniformDOS/'+filename+'_uniform_DOS_Zr.txt','w') as f:
    f.write('Spacing = '+str(uniform_X[1]-uniform_X[0])+'\n\n')
    f.write('Length  = '+str(len(uniform_X))+'\n\n')
    for i in range(0,len(uniform_Y)-6,6):
        for j in range(6):
            f.write(str('%.4E'%uniform_Y[i+j])+" ")
        f.write('\n')

"""




"""
num = 200
uniform_X = np.linspace(x[0],x[-1],num)
uniform_Y = [0.0]*num

for i in range(len(uniform_X)):
    if i > split:
        uniform_Y[i] = interpolate(x,y,uniform_X[i])
    else: 
        uniform_Y[i] = 0.0


inv_area = 1.0/np.trapz(uniform_Y,x=uniform_X)
uniform_Y = [val*inv_area for val in uniform_Y]

with open('uniformDOS/'+filename+'_uniform_DOS_H.txt','w') as f:
    f.write('Spacing = '+str(uniform_X[1]-uniform_X[0])+'\n\n')
    f.write('Length  = '+str(len(uniform_X))+'\n\n')
    for i in range(0,len(uniform_Y)-6,6):
        for j in range(6):
            f.write(str('%.4E'%uniform_Y[i+j])+" ")
        f.write('\n')








plt.plot(uniform_X,uniform_Y)
plt.show()



"""













