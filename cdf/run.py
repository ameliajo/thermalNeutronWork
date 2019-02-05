from test09_sab import *
import matplotlib.pyplot as plt
from plotHelp import *


def getSABval(sab,a,b,n_beta):
    return sab[a*n_beta+b]


def int_SAB_da(sab,b,alphas,n_beta):
    """
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5,0,0,0,0, 4,0,0,0,0, 1,0,0,0,0, 2,0,0,0,0]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=5)
    '0.085'

    >>> sab = [5,1,1,1,1, 4,1,1,1,1, 1,1,1,1,1, 2,1,1,1,1]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=5)
    '0.085'

    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> ['%.3f' % int_SAB_da(sab,b=x,alphas=alphas,n_beta=3) for x in range(3)]
    ['0.075', '0.069', '0.063']

    >>> sab = [5,1,1, 4,1,1, 1,1,1, 2,1,1]
    >>> alphas = [0.01,0.02,0.04,0.08]
    >>> '%.3f' % int_SAB_da(sab,b=0,alphas=alphas,n_beta=3)
    '0.155'
    """
    integral = 0.0
    for a in range(len(alphas)-1):
        yL = getSABval(sab,a,  b,n_beta)
        yR = getSABval(sab,a+1,b,n_beta)
        integral += (yL+yR)*0.5*(alphas[a+1]-alphas[a])
    return integral



def getEq18(sab,alphas,n_beta):
    """
    This solves Eq. 18 in Pavlou's paper  h*(a|b) = S(a,b) / int S(a,b) da
    So first we calculate the denominator, which is just a func of beta. Then
    we divide S(a,b) by that. The output is a list of lists, where eq18[a][b]
    gives h*(a|b)

    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> ['%.3f' % getEq18(sab,alphas,3)[0][b] for b in range(3)] # h*(a0|b0...2)
    ['66.667', '69.565', '73.016']
    >>> ['%.3f' % getEq18(sab,alphas,3)[1][b] for b in range(3)] # h*(a1|b0...2)
    ['13.333', '11.594', '9.524']
    >>> ['%.3f' % getEq18(sab,alphas,3)[2][b] for b in range(3)] # h*(a2|b0...2)
    ['26.667', '26.087', '25.397']
    >>> ['%.3f' % getEq18(sab,alphas,3)[3][b] for b in range(3)] # h*(a3|b0...2)
    ['53.333', '55.072', '57.143']
    """
    denominator = [int_SAB_da(sab,b,alphas,n_beta) for b in range(n_beta)]
    eq18 = [ [ getSABval(sab,a,b,n_beta)/denominator[b] \
                 for b in range(n_beta) ] \
                   for a in range(len(alphas)) ]
    return eq18


def generateEq19(sab,alphas,n_beta):
    """
    >>> alphas = [0.01,0.02,0.03,0.04]
    >>> sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
    >>> # H(a0,b0)=0   H(a0,b1)=0   H(a0,b2)=0 ---- CDF = 0 at first point
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[0]]
    ['0.000', '0.000', '0.000']
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[1]]
    ['0.400', '0.406', '0.413']
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[2]]
    ['0.600', '0.594', '0.587']
    >>> # H(a3,b0)=1   H(a3,b1)=1   H(a3,b2)=1 ---- CDF = 1 at final point
    >>> ['%.3f' % x for x in generateEq19(sab,alphas,n_beta=3)[3]]
    ['1.000', '1.000', '1.000']
    """
    eq19 = [[0.0]*n_beta for a in range(len(alphas))]
    for b in range(n_beta):
        eq18 = getEq18(sab,alphas,n_beta)
        for a in range(1,len(alphas)):
            delta = alphas[a]-alphas[a-1]
            eq19[a][b] = eq19[a-1][b] + (eq18[a-1][b]+eq18[a][b]) * 0.5 * delta
    return eq19

def generateEq20(sab,alphas,n_beta,a_min,a_max):
    eq19 = generateEq19(sab,alphas,n_beta)
    H  = [[0.0]*n_beta for a in range(len(alphas))]
    for b in range(len(betas)):
        H_min = eq19[a_min][b]
        H_max = eq19[a_max][b]
        for a in range(len(alphas)):
            H[a][b] = (eq19[a][b] - H_min) / (H_max - H_min)
    return H

"""
alphas = [0.01,0.02,0.03,0.04]
betas  = [0.1,0.2,0.3]
sab = [5.0,4.8,4.6,  1.0,0.8,0.6,   2.0,1.8,1.6,   4.0,3.8,3.6]
"""
n_beta = len(betas)
n_alpha = len(alphas)

 
H = generateEq20(sab,alphas,n_beta,0,n_alpha-1)
H_transpose = list(map(list, zip(*H)))

h = getEq18(sab,alphas,n_beta)
h_transpose = list(map(list, zip(*h)))

"""
for b in range(n_beta):
    plt.plot(alphas,H_transpose[b])
    plt.plot(alphas,H_transpose[b],'o')
"""

#for b in range(n_beta):
#    plt.plot(alphas,[getSABval(sab,a,b,n_beta) for a in range(n_alpha)])
#    plt.plot(alphas,[getSABval(sab,a,b,n_beta) for a in range(n_alpha)],'o')
#plt.show()


#scalarMap, colorBar = prepPlot(alphas)
#for b in range(n_beta):
#    plt.plot(alphas,H_transpose[b],label='beta: '+str(betas[b]),color=scalarMap.to_rgba(b))
#finishPlotting(colorBar,'H(a|b)')

scalarMap, colorBar = prepPlot(alphas)
for b in range(n_beta):
    plt.plot(alphas,h_transpose[b],label='beta: '+str(betas[b]),color=scalarMap.to_rgba(b))
finishPlotting(colorBar,'h(a|b)')





if __name__ == "__main__":
    import doctest
    doctest.testmod()
