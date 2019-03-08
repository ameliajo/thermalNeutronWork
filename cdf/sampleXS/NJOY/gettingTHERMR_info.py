import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def prepPlot(alphas):
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    colorBar = plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)
    plt.clf()
    return scalarMap, colorBar



def finishPlotting(colorBar):
    ax = plt.gca()
    plt.colorbar(colorBar).ax.set_ylabel('energies')
    plt.show()





def makeFloat(string):
    if string[-3] == '+': return float(string[:-3]+'E'+string[-3:]) 
    if string[-2] == '+': return float(string[:-2]+'E'+string[-2:])
    if string[-3] == '-': return float(string[:-3]+'E'+string[-3:]) 
    if string[-2] == '-': return float(string[:-2]+'E'+string[-2:])
    return float(string)



f = open('file6.txt')
line = f.readline()

allEnergies = []
allDatas = []
currentChunk = ''
startingNewChunk = True
n  = 0

while line:
    #print(line)
    vec = line.split()
    if len(vec) > 0:
        if vec[0] == '0.000000+0' and startingNewChunk:
            allEnergies.append(makeFloat(vec[1]))
            startingNewChunk = False
            if currentChunk != '':
                allDatas.append(currentChunk)
            currentChunk = ''
        else:
            currentChunk += ' '+line[:66]
            startingNewChunk = True
            

    line = f.readline()
    n += 1

allDatas.append(currentChunk)
f.close()

def makeVec(data):
    dataVec = data.split()
    finalDataVec = []
    for i in range(len(dataVec)-1):
        if len(dataVec[i]) >11:
            if dataVec[i][11] == '-':
                A = dataVec[i][:11]
                B = dataVec[i][11:]
                finalDataVec.append(makeFloat(A))
                finalDataVec.append(makeFloat(B))
            elif dataVec[i][10] == '-':
                A = dataVec[i][:10]
                B = dataVec[i][10:]
                finalDataVec.append(makeFloat(A))
                finalDataVec.append(makeFloat(B))
        else:
            finalDataVec.append(makeFloat(dataVec[i]))

    return finalDataVec




energies = allEnergies[:]
datas = allDatas[:]
scalarMap, colorBar = prepPlot(allEnergies)
for i,E_in in enumerate(energies):
    x,y,z = [],[],[]
    data = makeVec(datas[i])

    for j in range(int(len(data)/3)):
        x.append(data[3*j])
        y.append(data[3*j+1])
        z.append(data[3*j+2])


    plt.plot(x,y,label=str(E_in),color=scalarMap.to_rgba(i))
plt.xscale('log')
plt.yscale('log')

finishPlotting(colorBar)

#plt.legend(loc='best')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()








	
