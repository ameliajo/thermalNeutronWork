import subprocess
from makeTest09Rho import *

subprocess.run(['ls'])


alphaVals = [0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7 ]
betaVals = [7.6, 7.7, 7.8, 7.9, 8, 8.05, 8.1, 8.15, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7 ]
 
writeRho('inputVals.txt',2,continRho,alphaVals,betaVals)
subprocess.run(['ls'])




