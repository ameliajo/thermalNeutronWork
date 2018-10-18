import matplotlib.pyplot as plt

def plotErrorVsBeta(alpha,beta,y_delta,n_delta,colors,name):
    plt.figure(num=None, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    counter = 0
    for a in range(len(alpha)):
        error = []
        index = len(beta)*a
        for b in range(len(beta)):
            y,n = y_delta[index+b], n_delta[index+b]
            error.append(100.0*abs(y-n)/y)
            #error.append(abs(y-n))
            #error.append(abs(y))
        plt.plot(beta,error,label="alpha: "+str(alpha[a]),color=colors[counter])
        counter += 1

    plt.title("Error of Delta vs. Continuous Representation of Single Peak")
    plt.xlabel("Beta Values")
    plt.ylabel("Rel. Error (%)")
    plt.legend(loc='center left')

    plt.savefig(name, bbox_inches='tight')
    plt.show()



def plotErrorVsAlpha(alpha,beta,y_delta,n_delta,colors,name):
    plt.figure(num=None, figsize=(8, 6), dpi=100, facecolor='w', edgecolor='k')
    counter = 0
    for b in range(len(beta)):
        error = []
        for a in range(len(alpha)):
            index = b+a*len(beta)
            y,n = y_delta[index], n_delta[index]
            error.append(100.0*abs(y-n)/y)
            #if ((100.0*abs(y-n)/y) > 10): print("     ",beta[b],y,n)
            #error.append(abs(y-n))
            #error.append(abs(y))
        plt.plot(alpha,error,label="beta: "+str(beta[b]),color=colors[counter])
        counter += 1

    plt.title("Error of Delta vs. Continuous Representation of Single Peak")
    plt.xlabel("Alpha Values")
    plt.ylabel("Rel. Error (%)")
    plt.legend(loc='center left')

    plt.savefig(name, bbox_inches='tight')
    plt.show()



