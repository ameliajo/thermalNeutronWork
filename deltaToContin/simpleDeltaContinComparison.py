import matplotlib.pyplot as plt

no_delta   = [ 8.259370e-3, 8.345264e-3, 8.431157e-3, 1.628052e-2, 1.644894e-2, 1.661735e-2 ]

with_delta = [ 8.322095e-3, 8.408637e-3, 8.495179e-3, 1.640384e-2, 1.657352e-2, 1.674320e-2 ]

numA = 2
numB = 3

a = 0

for a in range(numA):
    vecToPlot_y_delta = []; vecToPlot_n_delta = []
    error = []
    for b in range(numB):
        n_delta = no_delta[numB*a+b]
        y_delta = with_delta[numB*a+b]
        error.append(abs(n_delta-y_delta)/y_delta)
    plt.plot(error)

plt.show()

# This is the one with no delta
#################################################################
# leapr
#  24 /
#  'delta functions without any contin'/
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
#  'delta functions without any contin'/
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
