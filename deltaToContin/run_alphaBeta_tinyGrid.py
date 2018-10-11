import matplotlib.pyplot as plt
from plotHelp import *

n_delta = [8.259370e-3,8.345264e-3,8.431157e-3,1.628052e-2,1.644894e-2,
           1.661735e-2]
y_delta = [8.322095e-3,8.408637e-3,8.495179e-3,1.640384e-2,1.657352e-2,
           1.674320e-2]
alpha   = [0.01,0.02]
beta    = [0.02,0.04,0.06]

colors  = [ "#ff0000", "#ffa600"]
plotErrorVsBeta(alpha,beta,y_delta,n_delta,colors,"alphaBeta_tinyGrid/plot_alphaBeta_tinyGrid.png")

colors  = [ "#ff0000","#ff5300","#ffa600"]
plotErrorVsAlpha(alpha,beta,y_delta,n_delta,colors,"alphaBeta_tinyGrid/plot_alphaBeta_tinyGrid.png")

# This is the one with no delta
#################################################################
# leapr
#  24 /
#  'simple peak with no delta functions, only contin'/
#  1 1/
#   101 1001/
#   0.99917 20.449 2 0 0/
#   1 1 15.85316 3.8883 1/
#  2 3 1/
#   0.01 0.02 /
#   0.02 0.04 0.06 /
#  296/
#   0.01 100/
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 101 1 1 1 1 1 1 1 1 1 1
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 /
#  0.0 0.0 10.0/  for diffusion, use 120. for second number
#   0/
#  ' h(h2o) thermal scattering '/
#  ' '/
#  ' temperatures = 296 deg k. '/
#  ' '/
#  ' shortened version of the endf/b-vi.4 evaluation for '/
#  ' hydrogen in water.  the energy transfer is limited to '/
#  ' 0.625 ev.  this is only for njoy testing, not for '/
#  ' real applications. '/
#  /
# stop

# This is the one with delta
#################################################################
# leapr
#  24 /
#  'simple peak with delta functions'/
#  1 1/
#   101 1001/
#   0.99917 20.449 2 0 0/
#   1 1 15.85316 3.8883 1/
#  2 3 1/
#   0.01 0.02 /
#   0.02 0.04 0.06 /
#  296/
#   0.01 100/
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
#   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 /
#  0.0 0.0 5.0/  for diffusion, use 120. for second number
#   1/
#   0.50/
#   5.0/
#  ' h(h2o) thermal scattering '/
#  ' '/
#  ' temperatures = 296 deg k. '/
#  ' '/
#  ' shortened version of the endf/b-vi.4 evaluation for '/
#  ' hydrogen in water.  the energy transfer is limited to '/
#  ' 0.625 ev.  this is only for njoy testing, not for '/
#  ' real applications. '/
#  /
# stop
