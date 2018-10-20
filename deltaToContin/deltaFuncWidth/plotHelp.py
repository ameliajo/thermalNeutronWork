import matplotlib.pyplot as plt

color1 = [ "#ff0000", "#ff0d00", "#ff1a00", "#ff2600", "#ff3300", "#ff4000", 
           "#ff4c00", "#ff5900", "#ff6600", "#ff7300", "#ff8000", "#ff8c00", 
           "#ff9900", "#ffa600", "#ffb200", "#ffbf00", "#ffcc00", "#ffd900", 
           "#ffe600", "#fff200", "#ffff00" ]

color2 = [ "#00cc00", "#00c20d", "#00b81a", "#00ad26", "#00a333", "#009940", 
           "#008f4c", "#008559", "#007a66", "#007073", "#006680", "#005c8c", 
           "#005299", "#0047a6", "#003db2", "#0033bf", "#0029cc", "#001fd9", 
           "#0014e6", "#000af2", "#0000ff" ]

color1shorter = [ "#ff0000", "#ff4c00", "#ff7300", "#ff9900", "#ffbf00", "#ffe600" ]
color2shorter = [ "#00cc00", "#008f4c", "#007073", "#005299", "#0033bf", "#0014e6" ]

def plotErrorAgainstBetaForGivenAlpha(a,yDelta,noDelta,alpha,beta):
    f = plt.figure(a)
    for i,n_delta in enumerate(noDelta):
        error = []
        #error2 = []
        index = len(beta)*a
        for b in range(len(beta)):
            y,n = yDelta[index+b], n_delta[index+b]
            #error.append(100.0*abs(y-n)/y)
            error.append(abs(y-n))
            #error.append(n)
            #if i == 0: error2.append(y)
        #if i == 0 :plt.plot(beta,error2,color=color2shorter[i+5],label="With Delta Function",linewidth=1)
        plt.plot(beta,error,label="Width: "+str(i+1)+" points used",color=color1shorter[i])

    plt.title(r"Run #2 NJOY-generated Error for $\alpha$ = "+str(alpha[a])+" comparing\ndelta input against various triangles to approximate delta function")
    #plt.ylabel("Relative Error in %")
    plt.ylabel("Absolute Error")
    plt.xlabel("Beta Values")
    #x1,x2,y1,y2 = plt.axis()
    #plt.axis((x1,x2,0.0,0.07))

    plt.legend(loc='upper right')

    f.show()





def plotSabAgainstBetaForGivenAlpha(a,yDelta,noDelta,alpha,beta):

    f2 = plt.figure(a)

    for i,ndelta in enumerate(noDelta):
        index = len(beta)*a
        deltaSab    = [yDelta[index+b] for b in range(len(beta))]
        triangleSab = [ndelta[index+b] for b in range(len(beta))]

        if i == 0: plt.plot(beta,deltaSab,color=color2shorter[i+5],label="With Delta Function",linewidth=1)
        plt.plot(beta,triangleSab,label="Width: "+str(i+1)+" points used",color=color1shorter[i])


    plt.title(r"Run #2 NJOY-generated $S(\alpha,\beta)$ for $\alpha$ = "+str(alpha[a])+" using conventional\ndelta input as well various triangles to approximate delta function")
    plt.ylabel("Freq. Dist. Values"); plt.xlabel("Beta Values")
    plt.legend(loc='center right')

    f2.show()




def plotJustWithDeltaFuncs(desiredAlpha,alpha,beta,yDelta):
    f1 = plt.figure(0)
    for a in [0,3,6,9,12]:
        index = int(len(beta))*a
        deltaSab = [yDelta[index+b] for b in range(len(beta))]
        plt.plot(beta,deltaSab,color=color2[a],label=r"$\alpha=$"+str(alpha[a]))
        plt.title("Run #1 NJOY-generated S(a,b) using conventional delta input")
        plt.ylabel("Freq. Dist. Values"); plt.xlabel("Beta Values")

        plt.legend(loc='center right')

    f1.show()





