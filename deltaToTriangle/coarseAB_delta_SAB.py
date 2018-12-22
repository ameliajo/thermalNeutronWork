import matplotlib.pyplot as plt 


nalpha = 6
nbeta  = 5

sab = [6.3582210947849514E-004, 6.6787473428553921E-004, 4.0170911240443454E-004, 1.4939387393537234E-004, 4.9973694134425670E-009, 6.2069547527636784E-003, 6.5193805752992067E-003, 3.9354522300335763E-003, 1.4951767270435831E-003, 7.3703965309386740E-007, 4.8823538479843720E-002, 5.1244887980880428E-002, 3.2070724539960141E-002, 1.4791001957050716E-002, 1.1541146037390215E-004, 8.5503457424853491E-002, 8.9509833521191795E-002, 6.5529907246419508E-002, 5.5178000555969269E-002, 3.1834459314046587E-003, 4.8322617692538039E-002, 5.0486021076118333E-002, 4.4349626444852366E-002, 5.9970597748557453E-002, 8.9613520067872918E-003, 8.0696197677929189E-005, 8.4739788279641258E-005, 1.2876224267362716E-004, 6.0204093788769402E-004, 9.5145843306852131E-003]



# Plotting all beta for given alpha
a_1_b_all = sab[0*nbeta:1*nbeta]
a_2_b_all = sab[1*nbeta:2*nbeta]
a_3_b_all = sab[2*nbeta:3*nbeta]
a_4_b_all = sab[3*nbeta:4*nbeta]
a_5_b_all = sab[4*nbeta:5*nbeta]
a_6_b_all = sab[5*nbeta:6*nbeta]

print(a_2_b_all)
betaVals = [0.0,0.1,1.0,5.0,25.0]
plt.plot(betaVals,a_1_b_all,label='alpha=0.01')
plt.plot(betaVals,a_2_b_all,label='alpha=0.10')
plt.plot(betaVals,a_3_b_all,label='alpha=1.00')
plt.plot(betaVals,a_4_b_all,label='alpha=5.00')
plt.plot(betaVals,a_5_b_all,label='alpha=10.0')
plt.plot(betaVals,a_6_b_all,label='alpha=50.0')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='best')
plt.show()





