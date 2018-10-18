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








def plotAgainstBetaForGivenAlpha(a,yDelta,noDelta,alpha,beta):
    f = plt.figure(a)
    for i,n_delta in enumerate(noDelta):
        error = []
        index = len(beta)*a
        for b in range(len(beta)):
            y,n = yDelta[index+b], n_delta[index+b]
            error.append(100.0*abs(y-n)/y)
        plt.plot(beta,error,label="Width: "+str(i+1)+" points used",color=color1shorter[i])

    plt.title("Error of Delta vs. Continuous Representation of Single Peak for Alpha = "+str(alpha[a]))
    plt.xlabel("Beta Values")
    plt.ylabel("Rel. Error (%)")
    plt.legend(loc='center right')

    #$plt.show()
    f.show()





