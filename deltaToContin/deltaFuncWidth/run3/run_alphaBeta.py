import matplotlib.pyplot as plt
from plotHelp import *
from run3_results import *

colors = [color1,color2]
alpha = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]
beta = [0.01,0.1,1.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0]


"""
yDelta_300K...yDelta_1300K
nDelta_300K...nDelta_1300K
"""

#plotErrorAgainstBetaForGivenAlpha(0, yDelta_300K,nDelta_300K,alpha,beta,color1shorter[0])
#plotErrorAgainstBetaForGivenAlpha(6, yDelta_300K,nDelta_300K,alpha,beta,color1shorter[2])
#plotErrorAgainstBetaForGivenAlpha(12,yDelta_300K,nDelta_300K,alpha,beta,color1shorter[4])
#for i in range(len(alpha)):
#    plotErrorAgainstBetaForGivenAlpha(i, yDelta_300K,nDelta_300K,alpha,beta,color1[i])

# alpha = 6.0
"""
plotErrorAgainstBetaForGivenAlpha(6, yDelta_300K,nDelta_300K,alpha,beta,color1shorter[0],300)
plotErrorAgainstBetaForGivenAlpha(6, yDelta_500K,nDelta_500K,alpha,beta,color1shorter[1],500)
plotErrorAgainstBetaForGivenAlpha(6, yDelta_700K,nDelta_700K,alpha,beta,color1shorter[2],700)
plotErrorAgainstBetaForGivenAlpha(6, yDelta_900K,nDelta_900K,alpha,beta,color1shorter[3],900)
plotErrorAgainstBetaForGivenAlpha(6, yDelta_1100K,nDelta_1100K,alpha,beta,color1shorter[4],1100)
plotErrorAgainstBetaForGivenAlpha(6, yDelta_1300K,nDelta_1300K,alpha,beta,color1shorter[5],1300)
"""
plotJustASab(12,yDelta_300K,nDelta_300K,color1shorter[0],300,alpha,beta)



plt.show()



#plotSabAgainstBetaForGivenAlpha(0,yDelta,noDelta,alpha,beta)
#plotSabAgainstBetaForGivenAlpha(6,yDelta,noDelta,alpha,beta)
#plotSabAgainstBetaForGivenAlpha(12,yDelta,noDelta,alpha,beta)

#plotErrorAgainstBetaForGivenAlpha(0,yDelta,noDelta,alpha,beta)
#plotErrorAgainstBetaForGivenAlpha(6,yDelta,noDelta,alpha,beta)
#plotErrorAgainstBetaForGivenAlpha(12,yDelta,noDelta,alpha,beta)

#desiredAlpha = [0,3,7,9,12]
#plotJustWithDeltaFuncs(desiredAlpha,alpha,beta,yDelta)

#input()
