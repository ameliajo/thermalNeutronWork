import matplotlib.pyplot as plt
from plotHelp import *
from run2_results import *

colors = [color1,color2]
alpha = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]
beta = [0.01,0.1,1.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0]


noDelta = [nDelta_width2,nDelta_width4,nDelta_width6,nDelta_width8,nDelta_width10]
#plotSabAgainstBetaForGivenAlpha(0,yDelta,noDelta,alpha,beta)
#plotSabAgainstBetaForGivenAlpha(6,yDelta,noDelta,alpha,beta)
#plotSabAgainstBetaForGivenAlpha(12,yDelta,noDelta,alpha,beta)

plotErrorAgainstBetaForGivenAlpha(0,yDelta,noDelta,alpha,beta)
plotErrorAgainstBetaForGivenAlpha(6,yDelta,noDelta,alpha,beta)
plotErrorAgainstBetaForGivenAlpha(12,yDelta,noDelta,alpha,beta)

#desiredAlpha = [0,3,7,9,12]
#plotJustWithDeltaFuncs(desiredAlpha,alpha,beta,yDelta)

input()
