import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def prepPlot(alphas):
    plt.clf()
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    return scalarMap, plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)


def makeFloat(string):
    if string[-3] == '+': return float(string[:-3]+'E'+string[-3:]) 
    if string[-2] == '+': return float(string[:-2]+'E'+string[-2:])
    if string[-3] == '-': return float(string[:-3]+'E'+string[-3:]) 
    if string[-2] == '-': return float(string[:-2]+'E'+string[-2:])
    return float(string)



def readFile3IntoVector():
    # This reads in from a txt document and separates file 6 into its differnt
    # initial energies
    f = open('file3.txt')
    line = f.readline()
    E_vec, xs_vec = [], []
    while line:
        vec = line.split()[:-2]
        if '1301' in vec[-1]:
            temp = vec[-1] + ' '
            vec[-1] = temp.replace('1301 ','')
        thisLine = [makeFloat(x) for x in vec]
        E_vec += [thisLine[i] for i in range(0,len(thisLine),2)]
        xs_vec += [thisLine[i] for i in range(1,len(thisLine),2)]
        line = f.readline()
    f.close()
    return E_vec, xs_vec 



def interpolate(xVec,yVec,x):
    if x < xVec[0]:  return yVec[0]
    if x > xVec[-1]: return yVec[-1]
    for i in range(len(xVec)-1):
        if xVec[i] <= x < xVec[i+1]:
            m = (yVec[i+1]-yVec[i])/(xVec[i+1]-xVec[i])
            b = yVec[i] - m*xVec[i]
            return m * x + b




def readFile6IntoVector():
    f = open('file6.txt')
    line = f.readline()
    allEnergies, allDatas = [], []
    currentChunk = ''
    startingNewChunk = True
    while line:
        vec = line.split()
        if len(vec) > 0:
            if vec[0] == '0.000000+0' and vec[2] == '0' and startingNewChunk:
                allEnergies.append(makeFloat(vec[1]))
                startingNewChunk = False
                if currentChunk != '': allDatas.append(currentChunk)
                currentChunk = ''
            else:
                currentChunk += ' '+line[:66]
                startingNewChunk = True
        line = f.readline()
    allDatas.append(currentChunk)
    f.close()
    return allEnergies,allDatas

def makeVecFile6(data):
    dataVec = data.split()
    finalDataVec = []
    for i in range(len(dataVec)):
        #print(dataVec[i])
        if len(dataVec[i]) > 11:
            minusLoc = None
            dotLoc   = None
            beginning = 0
            for j in range(len(dataVec[i])):
                #print('--------- ',dataVec[i][j])
                if dataVec[i][j] == '-': minusLoc = j
                if dataVec[i][j] == '.': dotLoc = j
                if minusLoc != None and dotLoc != None and dotLoc == minusLoc + 2:
                    #print(dataVec[i][beginning:minusLoc])
                    finalDataVec.append(makeFloat(dataVec[i][beginning:minusLoc]))
                    beginning = minusLoc
                    minusLoc = None
                    dotLoc = None
            finalDataVec.append(makeFloat(dataVec[i][beginning:]))
                
            #print(finalDataVec[-2:])
            #if dataVec[i][11] == '-':
            #    for j in range(0,len(dataVec[i])+1,11):
            #        if len(dataVec[i][j:j+11]) > 0:
            #            finalDataVec.append(makeFloat(dataVec[i][j:j+11]))
            #elif dataVec[i][10] == '-':
            #    for j in range(0,len(dataVec[i])+1,10):
            #        print(dataVec[i][j:j+10])
            #        if len(dataVec[i][j:j+10]) > 0:
            #            finalDataVec.append(makeFloat(dataVec[i][j:j+10]))
        else:
            finalDataVec.append(makeFloat(dataVec[i]))
    return finalDataVec




# ----------------------------------------------------------------------------
# Read in File3 XS Data
# ----------------------------------------------------------------------------
E_vec, xs_vec = readFile3IntoVector()
E_vec = E_vec[:-1]
xs_vec = xs_vec[:-1]
#plt.plot(E_vec,xs_vec)
#plt.xscale('log')
#plt.yscale('log')
#plt.show()


# ----------------------------------------------------------------------------
# Read in File6 E->E matrix
# ----------------------------------------------------------------------------
energies, datas = readFile6IntoVector()
#scalarMap, colorBar = prepPlot(energies)
#makeVecFile6(datas[1])
#for i,E_in in enumerate(energies):
#    x,y,z = [],[],[]
#    data = makeVecFile6(datas[i])
#    for j in range(int(len(data)/3)):
#        x.append(data[3*j])
#        y.append(data[3*j+1])
#        z.append(data[3*j+2])
#    plt.plot(x,y,label=str(E_in),color=scalarMap.to_rgba(i))
#plt.xscale('log')
#plt.yscale('log')
#ax = plt.gca()
#plt.colorbar(colorBar).ax.set_ylabel('energies')
#plt.show()


# ----------------------------------------------------------------------------
# Combine them both 
# ----------------------------------------------------------------------------
scalarMap, colorBar = prepPlot(energies)
for i,E_in in enumerate(energies):
    if E_in > 1e-4: continue
    x,y,z = [],[],[]
    file6_data = makeVecFile6(datas[i])
    xsValue = interpolate(E_vec,xs_vec,E_in)
    for j in range(int(len(file6_data)/3)):
        x.append(file6_data[3*j])
        y.append(file6_data[3*j+1]*xsValue)
    plt.plot(x,y,label=str(E_in),color=scalarMap.to_rgba(i))
plt.xscale('log')
plt.yscale('log')
ax = plt.gca()
plt.colorbar(colorBar).ax.set_ylabel('energies')
plt.show()


"""


"""
	
