from numpy import exp
import numpy as np
import matplotlib.pyplot as plt



def interpolate(xList,yList,x):
    for i in range(len(xList)-1):
        if (xList[i] <= x < xList[i+1]):
            m = (yList[i+1]-yList[i]) / (xList[i+1]-xList[i])
            b = yList[i] - m*xList[i]
            return m*x+b
    return 0.0 if (xList[-1] < x) or xList[0] > x else yList[-1]


inputBetas = list(np.linspace(0.0,4,10)) + list(np.linspace(1.0,3,25))
inputBetas.sort()
inputT1 = [-2*x*x+8*x for x in inputBetas]



# REFLECT BETA AND T1
reflectT1 = [interpolate(inputBetas,inputT1,abs(beta)) for beta in reflectBetas]
reflectBetas = np.linspace(-10,10,21) #reflectBetas = np.linspace(-2,8,11)
plt.plot(inputBetas,inputT1,'r')
plt.plot(inputBetas,inputT1,'ro',markersize=5)
plt.plot(reflectBetas,reflectT1,'b')
plt.plot(reflectBetas,reflectT1,'bo',markersize=1)

for i in range(int(len(reflectBetas)/2)+1):
    reflectT1[i] *= exp(-reflectBetas[i])

plt.plot(reflectBetas,reflectT1,'g')
plt.plot(reflectBetas,reflectT1,'go',markersize=2)
plt.show()
"""
"""

#betas = np.linspace(0,10,21)
betas = [x for x in reflectBetas if x >= 0]
T1 = [interpolate(inputBetas,inputT1,beta) for beta in betas]

T2 = [0.0]*len(betas)

for b in range(len(betas)):
    beta = betas[b]
    for b_prime in range(len(reflectBetas[:-1])):
        betaPrime = reflectBetas[b_prime]
        L = interpolate(betas,T1,abs(betaPrime)) * \
            interpolate(betas,T1,abs(beta-reflectBetas[b_prime]))
        R = interpolate(betas,T1,abs(reflectBetas[b_prime+1])) * \
            interpolate(betas,T1,abs(beta-reflectBetas[b_prime+1]))
        if      reflectBetas[b_prime]   < 0: L*= exp(-reflectBetas[b_prime])
        if beta-reflectBetas[b_prime]   < 0: L*= exp(-beta+reflectBetas[b_prime])
        if      reflectBetas[b_prime+1] < 0: R*= exp(-reflectBetas[b_prime+1])
        if beta-reflectBetas[b_prime+1] < 0: R*= exp(-beta+reflectBetas[b_prime+1])

        T2[b] += (L+R)*(reflectBetas[b_prime+1]-reflectBetas[b_prime])*0.5

plt.plot(betas,T2,'k')
plt.plot(betas,T2,'ko')
#plt.show()




"""







plt.show()



#move T1 to ultrafine grid
T1 = [0.0]*len(reflectBetas)
for b,beta in enumerate(reflectBetas):
    T1[b] = interpolate(inputBetas,inputT1,beta)

plt.plot(inputBetas,inputT1,'r')
plt.plot(inputBetas,inputT1,'ro',markersize=5)
plt.plot(reflectBetas,T1,'b')
plt.plot(reflectBetas,T1,'bo',markersize=1)


T2 = [0.0]*len(reflectBetas)

for b,beta in enumerate(reflectBetas):
    for b_prime,beta_prime in enumerate(reflectBetas[:-1]):
        L = interpolate(reflectBetas,T1,reflectBetas[b_prime]) * \
            interpolate(reflectBetas,T1,beta-reflectBetas[b_prime])
        R = interpolate(reflectBetas,T1,reflectBetas[b_prime+1]) * \
            interpolate(reflectBetas,T1,beta-reflectBetas[b_prime+1])

        print(beta,beta_prime,L,R)
        T2[b] += (L+R)*(reflectBetas[b_prime+1]-reflectBetas[b_prime])*0.5

plt.plot(reflectBetas,T2,'g')
plt.plot(reflectBetas,T2,'go',markersize=1)
plt.show()
"""








#betas = [-2.0, 0.0, 2.0]
#T1 = [3.0*exp(2),1.0,3.0]



