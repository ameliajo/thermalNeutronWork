import matplotlib.pyplot as plt
from plotHelp import *
from run1_results import *

color1 = [ "#ff0000", "#ff0d00", "#ff1a00", "#ff2600", "#ff3300", "#ff4000", 
           "#ff4c00", "#ff5900", "#ff6600", "#ff7300", "#ff8000", "#ff8c00", 
           "#ff9900", "#ffa600", "#ffb200", "#ffbf00", "#ffcc00", "#ffd900", 
           "#ffe600", "#fff200", "#ffff00" ]
color2 = [ "#00cc00", "#00c20d", "#00b81a", "#00ad26", "#00a333", "#009940", 
           "#008f4c", "#008559", "#007a66", "#007073", "#006680", "#005c8c", 
           "#005299", "#0047a6", "#003db2", "#0033bf", "#0029cc", "#001fd9", 
           "#0014e6", "#000af2", "#0000ff" ]
colors = [color1,color2]
alpha = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]
beta = [0.01,0.1,1.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0]

noDelta = [nDelta_width2,nDelta_width4,nDelta_width6,nDelta_width8,nDelta_width10]
plotAgainstBetaForGivenAlpha(0,yDelta,noDelta,alpha,beta)
plotAgainstBetaForGivenAlpha(6,yDelta,noDelta,alpha,beta)
plotAgainstBetaForGivenAlpha(12,yDelta,noDelta,alpha,beta)
input()
