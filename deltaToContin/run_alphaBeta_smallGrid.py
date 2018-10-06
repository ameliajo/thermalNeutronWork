import matplotlib.pyplot as plt
from plotHelp import *

colors = [ "#ff0000", "#ff0d00", "#ff1a00", "#ff2600", "#ff3300", "#ff4000", 
           "#ff4c00", "#ff5900", "#ff6600", "#ff7300", "#ff8000", "#ff8c00", 
           "#ff9900", "#ffa600", "#ffb200", "#ffbf00", "#ffcc00", "#ffd900", 
           "#ffe600", "#fff200", "#ffff00" ]
n_delta = [ 8.259370e-3, 8.345264e-3, 8.431157e-3, 8.517050e-3, 8.602944e-3, 
  8.688837e-3, 8.774731e-3, 8.860624e-3, 8.946517e-3, 9.032411e-3, 9.118304e-3, 
  9.204198e-3, 9.290091e-3, 9.375984e-3, 9.461878e-3, 1.628052e-2, 1.644894e-2, 
  1.661735e-2, 1.678577e-2, 1.695419e-2, 1.712261e-2, 1.729103e-2, 1.745945e-2, 
  1.762786e-2, 1.779628e-2, 1.796470e-2, 1.813312e-2, 1.830154e-2, 1.846996e-2, 
  1.863838e-2, 2.406878e-2, 2.431646e-2, 2.456414e-2, 2.481182e-2, 2.505950e-2, 
  2.530717e-2, 2.555485e-2, 2.580253e-2, 2.605021e-2, 2.629789e-2, 2.654557e-2, 
  2.679325e-2, 2.704092e-2, 2.728860e-2, 2.753628e-2, 3.162938e-2, 3.195315e-2, 
  3.227693e-2, 3.260070e-2, 3.292448e-2, 3.324825e-2, 3.357203e-2, 3.389580e-2, 
  3.421958e-2, 3.454335e-2, 3.486713e-2, 3.519090e-2, 3.551468e-2, 3.583845e-2, 
  3.616222e-2, 3.896742e-2, 3.936423e-2, 3.976103e-2, 4.015784e-2, 4.055464e-2, 
  4.095145e-2, 4.134825e-2, 4.174506e-2, 4.214186e-2, 4.253867e-2, 4.293547e-2, 
  4.333228e-2, 4.372908e-2, 4.412589e-2, 4.452269e-2, 4.608793e-2, 4.655479e-2, 
  4.702165e-2, 4.748852e-2, 4.795538e-2, 4.842225e-2, 4.888911e-2, 4.935597e-2, 
  4.982284e-2, 5.028970e-2, 5.075657e-2, 5.122343e-2, 5.169029e-2, 5.215716e-2, 
  5.262402e-2, 5.299580e-2, 5.352984e-2, 5.406389e-2, 5.459793e-2, 5.513198e-2, 
  5.566602e-2, 5.620007e-2, 5.673411e-2, 5.726816e-2, 5.780220e-2, 5.833625e-2, 
  5.887029e-2, 5.940434e-2, 5.993838e-2, 6.047243e-2, 5.969584e-2, 6.029428e-2, 
  6.089272e-2, 6.149115e-2, 6.208959e-2, 6.268803e-2, 6.328647e-2, 6.388490e-2, 
  6.448334e-2, 6.508178e-2, 6.568022e-2, 6.627865e-2, 6.687709e-2, 6.747553e-2, 
  6.807397e-2, 6.619277e-2, 6.685290e-2, 6.751303e-2, 6.817316e-2, 6.883329e-2, 
  6.949342e-2, 7.015354e-2, 7.081367e-2, 7.147380e-2, 7.213393e-2, 7.279406e-2, 
  7.345419e-2, 7.411432e-2, 7.477445e-2, 7.543458e-2, 7.249119e-2, 7.321040e-2, 
  7.392960e-2, 7.464881e-2, 7.536802e-2, 7.608722e-2, 7.680643e-2, 7.752563e-2, 
  7.824484e-2, 7.896404e-2, 7.968325e-2, 8.040245e-2, 8.112166e-2, 8.184086e-2, 
  8.256007e-2 ]

y_delta = [ 8.322095e-3, 8.408637e-3, 8.495179e-3, 8.581722e-3, 8.668264e-3, 
  8.754806e-3, 8.841348e-3, 8.927890e-3, 9.014433e-3, 9.100975e-3, 9.187517e-3, 
  9.274059e-3, 9.360602e-3, 9.447144e-3, 9.533686e-3, 1.640384e-2, 1.657352e-2, 
  1.674320e-2, 1.691288e-2, 1.708256e-2, 1.725224e-2, 1.742192e-2, 1.759160e-2, 
  1.776128e-2, 1.793096e-2, 1.810064e-2, 1.827032e-2, 1.844000e-2, 1.860968e-2, 
  1.877936e-2, 2.425062e-2, 2.450014e-2, 2.474966e-2, 2.499918e-2, 2.524870e-2, 
  2.549822e-2, 2.574774e-2, 2.599726e-2, 2.624678e-2, 2.649630e-2, 2.674582e-2, 
  2.699534e-2, 2.724486e-2, 2.749438e-2, 2.774390e-2, 3.186772e-2, 3.219389e-2, 
  3.252005e-2, 3.284621e-2, 3.317237e-2, 3.349854e-2, 3.382470e-2, 3.415086e-2, 
  3.447702e-2, 3.480319e-2, 3.512935e-2, 3.545551e-2, 3.578167e-2, 3.610783e-2, 
  3.643400e-2, 3.926031e-2, 3.966002e-2, 4.005972e-2, 4.045943e-2, 4.085914e-2, 
  4.125885e-2, 4.165855e-2, 4.205826e-2, 4.245797e-2, 4.285768e-2, 4.325738e-2, 
  4.365709e-2, 4.405680e-2, 4.445651e-2, 4.485621e-2, 4.643344e-2, 4.690370e-2, 
  4.737395e-2, 4.784420e-2, 4.831445e-2, 4.878470e-2, 4.925495e-2, 4.972521e-2, 
  5.019546e-2, 5.066571e-2, 5.113596e-2, 5.160621e-2, 5.207647e-2, 5.254672e-2, 
  5.301697e-2, 5.339209e-2, 5.392998e-2, 5.446786e-2, 5.500575e-2, 5.554364e-2, 
  5.608153e-2, 5.661942e-2, 5.715731e-2, 5.769520e-2, 5.823309e-2, 5.877098e-2, 
  5.930886e-2, 5.984675e-2, 6.038464e-2, 6.092253e-2, 6.014110e-2, 6.074381e-2, 
  6.134652e-2, 6.194923e-2, 6.255194e-2, 6.315465e-2, 6.375736e-2, 6.436007e-2, 
  6.496278e-2, 6.556549e-2, 6.616820e-2, 6.677091e-2, 6.737362e-2, 6.797633e-2, 
  6.857904e-2, 6.668523e-2, 6.735004e-2, 6.801484e-2, 6.867964e-2, 6.934445e-2, 
  7.000925e-2, 7.067406e-2, 7.133886e-2, 7.200366e-2, 7.266847e-2, 7.333327e-2, 
  7.399808e-2, 7.466288e-2, 7.532769e-2, 7.599249e-2, 7.302915e-2, 7.375341e-2, 
  7.447767e-2, 7.520192e-2, 7.592618e-2, 7.665044e-2, 7.737470e-2, 7.809895e-2, 
  7.882321e-2, 7.954747e-2, 8.027172e-2, 8.099598e-2, 8.172024e-2, 8.244449e-2, 
  8.316875e-2 ]

alpha = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
beta  = [0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30]

plotErrorVsBeta(alpha,beta,y_delta,n_delta,colors,"alphaBeta_smallGrid/plot_alphaBeta_smallGrid.png")
plotErrorVsAlpha(alpha,beta,y_delta,n_delta,colors,"alphaBeta_smallGrid/plot_alphaBeta_smallGrid.png")

# This is the one with no delta
#################################################################
#leapr
# 24 /
# 'simple peak with no delta functions, only contin'/
# 1 1/
#  101 1001/
#  0.99917 20.449 2 0 0/
#  1 1 15.85316 3.8883 1/
# 10 15 1/
#  0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10/
#  0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20
#  0.22 0.24 0.26 0.28 0.30/
# 296/
#  0.01 100/
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 101 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 /
# 0.0 0.0 10.0/  for diffusion, use 120. for second number
#  0/
# ' h(h2o) thermal scattering '/
# ' '/
# ' temperatures = 296 deg k. '/
# ' '/
# ' shortened version of the endf/b-vi.4 evaluation for '/
# ' hydrogen in water.  the energy transfer is limited to '/
# ' 0.625 ev.  this is only for njoy testing, not for '/
# ' real applications. '/
# /
#stop


# This is the one with delta
#################################################################
#leapr
# 24 /
# 'simple peak with delta functions'/
# 1 1/
#  101 1001/
#  0.99917 20.449 2 0 0/
#  1 1 15.85316 3.8883 1/
# 10 15 1/
#  0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10/
#  0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20
#  0.22 0.24 0.26 0.28 0.30/
# 296/
#  0.01 100/
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 /
# 0.0 0.0 5.0/  for diffusion, use 120. for second number
#  1/
#  0.50/
#  5.0/
# ' h(h2o) thermal scattering '/
# ' '/
# ' temperatures = 296 deg k. '/
# ' '/
# ' shortened version of the endf/b-vi.4 evaluation for '/
# ' hydrogen in water.  the energy transfer is limited to '/
# ' 0.625 ev.  this is only for njoy testing, not for '/
# ' real applications. '/
# /
#stop
