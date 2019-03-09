import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np


def prepPlot(alphas):
    plt.clf()
    cnorm = colors.Normalize(vmin=0,vmax=len(alphas)+2)
    scalarMap = cmx.ScalarMappable(norm=cnorm,cmap=plt.get_cmap('tab20')) #hot autumn tab10
    mymap = colors.LinearSegmentedColormap.from_list('funTestColors',\
            [scalarMap.to_rgba(a) for a in range(len(alphas))])
    return scalarMap, plt.contourf([[0,0],[0,0]], alphas, cmap=mymap)

def makeFloat(string):
    if string[-3] in '+-': return float(string[:-3]+'E'+string[-3:]) 
    if string[-2] in '+-': return float(string[:-2]+'E'+string[-2:])
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
            vec[-1] = str(vec[-1]+' ').replace('1301 ','')
        E_vec  += [makeFloat(vec[i]) for i in range(0,len(vec),2)]
        xs_vec += [makeFloat(vec[i]) for i in range(1,len(vec),2)]
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
    thisChunk = ''
    newChunk = True
    while line:
        vec = line.split()
        if len(vec) > 0:
            if vec[0] == '0.000000+0' and vec[2] == '0' and newChunk:
                allEnergies.append(makeFloat(vec[1]))
                if thisChunk != '': allDatas.append(thisChunk)
                thisChunk = ''
                newChunk = False
            else:
                thisChunk += ' '+line[:66]
                newChunk = True
        line = f.readline()
    f.close()
    allDatas.append(thisChunk)
    return allEnergies,allDatas

def makeVecFile6(data):
    dataVec = data.split()
    finalDataVec = []
    for i in range(len(dataVec)):
        beginning = 0
        if len(dataVec[i]) > 11:
            # These pieces are to account for negative numbers that wont work 
            # with x.split()
            minusLoc, dotLoc = None, None
            for j in range(len(dataVec[i])):
                if dataVec[i][j] == '-': minusLoc = j
                if dataVec[i][j] == '.': dotLoc = j
                if minusLoc != None and dotLoc != None and dotLoc == minusLoc + 2:
                    finalDataVec.append(makeFloat(dataVec[i][beginning:minusLoc]))
                    beginning = minusLoc
                    minusLoc, dotLoc = None, None
        finalDataVec.append(makeFloat(dataVec[i][beginning:]))
    return finalDataVec










if __name__ == '__main__':

    toPrint = {'xs':True, 'probMatrix':True, 'fullScat':True}
    toPrint = {'xs':False, 'probMatrix':False, 'fullScat':False, 'specificE':True}

    # --------------------------------------------------------------------------
    # Read in File3 XS Data
    # --------------------------------------------------------------------------
    E_vec, xs_vec = readFile3IntoVector()
    E_vec = E_vec[:-1]
    xs_vec = xs_vec[:-1]


    # --------------------------------------------------------------------------
    # Read in File6 P(E->E') matrix
    # --------------------------------------------------------------------------
    energies, datas = readFile6IntoVector()
    scalarMap, colorBar = prepPlot(energies)
    makeVecFile6(datas[1])


    if toPrint['xs']:
        plt.plot(E_vec,xs_vec)
        plt.xscale('log'); plt.yscale('log')
        plt.show()

    if toPrint['probMatrix']:
        for i,E_in in enumerate(energies):
            x,y = [],[]
            data = makeVecFile6(datas[i])
            for j in range(int(len(data)/3)):
                x.append(data[3*j])
                y.append(data[3*j+1])
            plt.plot(x,y,label=str(E_in),color=scalarMap.to_rgba(i))
        plt.xscale('log'); plt.yscale('log')
        plt.colorbar(colorBar).ax.set_ylabel('energies')
        ax = plt.gca()
        plt.show()


    if toPrint['fullScat']:
        areas = []
        scalarMap, colorBar = prepPlot(energies)
        for i,E_in in enumerate(energies):
            x,y = [],[]
            file6_data = makeVecFile6(datas[i])
            xsValue = interpolate(E_vec,xs_vec,E_in)
            for j in range(int(len(file6_data)/3)):
                x.append(file6_data[3*j])
                y.append(file6_data[3*j+1]*xsValue)
            plt.plot(x,y,label=str(E_in),color=scalarMap.to_rgba(i))
        plt.xscale('log'); plt.yscale('log')
        plt.colorbar(colorBar).ax.set_ylabel('energies')
        ax = plt.gca()
        plt.show()

    if toPrint['specificE']:

        specifiedE = [0.0005,0.0253,0.2907]

        for E in specifiedE:
            index = None
            for i in range(len(energies)-1):
                if energies[i] <= E < energies[i+1]:
                    index = i; break
            if index != None:
                x,y = [],[]
                file6_data = makeVecFile6(datas[index])
                xsValue = interpolate(E_vec,xs_vec,energies[index])
                for j in range(int(len(file6_data)/3)):
                    x.append(file6_data[3*j])
                    y.append(file6_data[3*j+1]*xsValue)

                print(energies[index])
                
                plt.plot(x,y,label=str(E)+' eV')
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel("E' [eV]")
        plt.ylabel("Prob / eV")
        plt.legend(loc='best')
        ax = plt.gca()
        plt.show()



