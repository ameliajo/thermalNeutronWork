import matplotlib.pyplot as plt 



def isValidABCombo(alpha,beta,A0,E,kbT):
    a_min = (E**0.5 - (E+beta*kbT)**0.5)**2 / (A0*kbT)
    a_max = (E**0.5 + (E+beta*kbT)**0.5)**2 / (A0*kbT)
    return a_min < alpha < a_max


def plotBetaForVariousAlpha(alphaVals,betaVals,sab,A0,E,kbT,scalarMap,style,addLabel):

    nalpha, nbeta = len(alphaVals), len(betaVals)
    # Plotting all beta for given alpha
    for a,alpha in enumerate(alphaVals):
        a_i_b_all = sab[a*nbeta:(a+1)*nbeta]
        validBeta = []
        validSab  = []
        for b,beta in enumerate(betaVals):
            if isValidABCombo(alpha,beta,A0,E,kbT):
                validBeta.append(beta)
                validSab.append(a_i_b_all[b])
        if addLabel:
            plt.plot(validBeta,validSab,label='alpha: '+str(alpha),color=scalarMap.to_rgba(a),marker=style)
        else:
            plt.plot(validBeta,validSab,color=scalarMap.to_rgba(a),marker=style,linestyle='--')


    #plt.legend(loc='best')
    plt.xlabel('beta')
    plt.ylabel('S(a,b)')


def plotErrorBetaForVariousAlpha(alphaVals,betaVals,sabGood,sabTest,A0,E,kbT,scalarMap):

    nalpha, nbeta = len(alphaVals), len(betaVals)
    # Plotting all beta for given alpha
    for a,alpha in enumerate(alphaVals):
        a_i_b_all_Good = sabGood[a*nbeta:(a+1)*nbeta]
        a_i_b_all_Test = sabTest[a*nbeta:(a+1)*nbeta]
        validBeta = []
        validSab  = []
        for b,beta in enumerate(betaVals):
            if isValidABCombo(alpha,beta,A0,E,kbT):
                validBeta.append(beta)
                validSab.append(100.0*(a_i_b_all_Test[b]-a_i_b_all_Good[b])/a_i_b_all_Good[b])
        plt.plot(validBeta,validSab,label='alpha: '+str(alpha),color=scalarMap.to_rgba(a))


    #plt.legend(loc='best')
    plt.xlabel('beta')
    plt.ylabel('S(a,b) Error (in %)')



def plotAlphaForVariousBeta(alphaVals,betaVals,sab,A0,E,kbT,scalarMap,style):

    nalpha, nbeta = len(alphaVals), len(betaVals)
    # Plotting all beta for given alpha
    for b,beta in enumerate(betaVals):
        b_i_a_all = [sab[b+(a*nbeta)] for a in range(nalpha)]
        validAlpha = []
        validSab   = []
        for a,alpha in enumerate(alphaVals):
            if isValidABCombo(alpha,beta,A0,E,kbT):
                validAlpha.append(alpha)
                validSab.append(b_i_a_all[a])
        plt.plot(validAlpha,validSab,label='beta: '+str(beta),color=scalarMap.to_rgba(b),linestyle=style)


    #plt.legend(loc='best')
    plt.xlabel('alpha')
    plt.ylabel('S(a,b)')






