#from test09_sab import *
import matplotlib.pyplot as plt


def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]


def int_SAB_da(sab,b,alphas,n_beta):
    """
    >>> sab = [5,0,0,0,0, 4,0,0,0,0, 1,0,0,0,0, 2,0,0,0,0]
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=5)
    '0.085'
    >>> sab = [5,1,1,1,1, 4,1,1,1,1, 1,1,1,1,1, 2,1,1,1,1]
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=5)
    '0.085'
    >>> sab = [5,1,1, 4,1,1, 1,1,1, 2,1,1]
    >>> alphas = [0.01,0.02,0.04,0.08]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=3)
    '0.155'
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=3)
    '0.075'
    >>> '%.3f' % int_SAB_da(sab,b=1,alphas=alphas,n_beta=3)
    '0.069'
    >>> '%.3f' % int_SAB_da(sab,b=2,alphas=alphas,n_beta=3)
    '0.063'
    """
    integral = 0.0
    for a in range(len(alphas)-1):
        yL = getSABval(sab,a,  b,n_beta)
        yR = getSABval(sab,a+1,b,n_beta)
        integral += (yL+yR)*0.5*(alphas[a+1]-alphas[a])
    return integral



def getEq18(sab,alphas,n_beta):
    """
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> ['%.3f' % getEq18(sab,alphas,3)[0][b] for b in range(3)]
    ['66.667', '69.565', '73.016']
    >>> ['%.3f' % getEq18(sab,alphas,3)[1][b] for b in range(3)]
    ['13.333', '11.594', '9.524']
    >>> ['%.3f' % getEq18(sab,alphas,3)[2][b] for b in range(3)]
    ['26.667', '26.087', '25.397']
    >>> ['%.3f' % getEq18(sab,alphas,3)[3][b] for b in range(3)]
    ['53.333', '55.072', '57.143']
    """
    denominator = [0.0]*n_beta
    for b in range(n_beta):
        denominator[b] = int_SAB_da(sab,b,alphas,n_beta)

    eq18 = [[0.0]*n_beta for a in range(len(alphas))]
    for a in range(len(alphas)):
        for b in range(n_beta):
            eq18[a][b] = getSABval(sab,a,b,n_beta)/denominator[b]

    return eq18

def generateEq19(sab,alphas,n_beta):
    """
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[0]]
    ['0.000', '0.000', '0.000']
    >>> '%.3f' % generateEq19(sab,alphas,n_beta=3)[1][0]
    '0.400'
    >>> '%.3f' % generateEq19(sab,alphas,n_beta=3)[2][0]
    '0.600'
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[3]]
    ['1.000', '1.000', '1.000']
    """
    #>>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[0][0]]
    #['0.000', '0.400', '0.600', '1.000']

    eq19 = [[0.0]*n_beta for a in range(len(alphas))]
    for b in range(n_beta):
        eq18 = getEq18(sab,alphas,n_beta)
        earlier = eq19[0][b]
        for a in range(1,n_alpha):
            delta = alphas[a]-alphas[a-1]
            eq19[a][b] = earlier + (eq18[a-1][b]+eq18[a][b]) * 0.5 * delta
            earlier = eq19[a][b]
    return eq19


alphas = [0.01,0.02,0.03,0.04]
betas  = [0.1,0.2,0.3]
sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]

betas  = [0.1]
sab = [5.0,1.0,2.0,4.0]
n_beta = len(betas)
n_alpha = len(alphas)


 

for b in range(n_beta):
    plt.plot(alphas,[getSABval(sab,a,b,n_beta) for a in range(n_alpha)])
    plt.plot(alphas,[getSABval(sab,a,b,n_beta) for a in range(n_alpha)],'o')






#plt.show()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
