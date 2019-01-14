from makeTest09Rho import * 
import subprocess



intro = \
    "leapr"+str('\n')+\
    "24 /"+str('\n')+\
    "'h in h2o, shortened endf model'/"+str('\n')+\
    "1 2/"+str('\n')+\
    " 101 1001/"+str('\n')+\
    " 0.99917 20.449 2 0 0/"+str('\n')+\
    " 1 1 15.85316 3.8883 1/"+str('\n')


def writeTillEndOfLine(f,vec):
    lineLen = 0
    for x in vec:
        if len(str(x))+lineLen < 70:
            f.write(str(x)+" ")
            lineLen += len(str(x))+1
        else:
            f.write("\n"+str(x)+" ")
            lineLen = len(str(x))
    f.write("/\n")


def generateNjoyInput(fileName,alphaVals,betaVals,phononDist,deltaFuncs):
    with open(fileName,'w') as f:
        f.write(intro)
        f.write(str(len(alphaVals))+" "+str(len(betaVals))+"/\n")
        writeTillEndOfLine(f,alphaVals)
        writeTillEndOfLine(f,betaVals)

        f.write("296/\n")
        f.write("0.00255 "+str(len(phononDist))+" /\n")
        writeTillEndOfLine(f,phononDist)

        if deltaFuncs:
            f.write(\
             "0.0 0. 0.5 / "+str('\n')+\
             " 2/"+str('\n')+\
             " .204 .4794/"+str('\n')+\
             " .166667 .333333/"+str('\n')+\
             "/"+str('\n')+\
            "stop"+str('\n'))
        else:
            f.write(\
             "0.0 0. 1.0 / "+str('\n')+\
             " 0/"+str('\n')+\
             "/"+str('\n')+\
            "stop"+str('\n'))


def runNJOY(fileName):
    subprocess.run(['cp',fileName,'/Users/amelia/NJOY2016/bin'])
    subprocess.call("~/NJOY2016/bin/njoy < "+str(fileName),shell=True)
    subprocess.run(['mv','./sab.txt','./NJOY_LEAPR/sabResults/sab_'+fileName+'.txt'])
    subprocess.run(['mv',fileName,'./NJOY_LEAPR/njoyInputs/'])
    #subprocess.run(['rm','tape24','output'])


if __name__=='__main__':

    alphaVals = [0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
    betaVals = [7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7 ]
 
    continRho = [0, .0005, .001, .002, .0035, .005, .0075, .01, .013, .0165, .02,  \
      .0245, .029, .034, .0395, .045, .0506, .0562, .0622, .0686, .075, .083, .091,\
      .099, .107, .115, .1197, .1214, .1218, .1195, .1125, .1065, .1005, .09542,   \
      .09126, .0871, .0839, .0807, .07798, .07574, .0735, .07162, .06974, .06804,  \
      .06652, .065, .0634, .0618, .06022, .05866, .0571, .05586, .05462, .0535,    \
      .0525, .0515, .05042, .04934, .04822, .04706, .0459, .04478, .04366, .04288, \
      .04244, .042, 0.0]

    width = 2
    fileName = 'triangleOfWidth'+str(width)
    generateNjoyInput(fileName,alphaVals,betaVals,getPhononDist(width,continRho),False)
    runNJOY(fileName)




