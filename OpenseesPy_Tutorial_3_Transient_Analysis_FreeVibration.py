# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 06:37:08 2023

@author: carlo
"""

import openseespy.opensees as op

op.wipe()
# General model definition (2 dimensions and 3 deegrees of freedom)
op.model('Basic', '-ndm', 2, '-ndf', 3)

########################### PHYSICAL PROPERTIES ###############################
# Shaft length
L = 1.5                  # (m)
# Shaft diameter
D = 0.5 * 0.0254         # (m)
# Lumped mass at 0.7L from the left end
m = 40                   # (kg)
import math
pi = math.acos(-1.0)
# Young´s Modulus
E = 2.1e11                # (Pa)
# Shaft Density
ro = 7850                 # (kg/m^3)
# Shaft cross-section area
A = 0.25 * pi * D**2      # (m^2)
# Shaft cross-section second moment of Inertia
Iz = pi * D**4 / 64       # (m^4)

############################### DEFINE NODES ##################################
# Total number of nodes
Nmax = 11
for i in range(1, Nmax + 1):
    op.node(i, (i - 1) * L / (Nmax - 1), 0)

############################## DOF CONSTRAINTS ################################
op.fix(   1, 1, 1, 0)
op.fix(Nmax, 1, 1, 0)

##################### ADDITION OF CONCENTRATION MASS ##########################
node_loaded = 8
op.mass(node_loaded, 0.0, m, 0.0)

############################# DEFINE ELEMENTS #################################
# Geometric transformation Tag
transTag = 1
op.geomTransf('Linear', transTag)
for index in range(1, Nmax):
    op.element('elasticBeamColumn', index, *[index, index + 1], A, E, Iz, 
               transTag, '-mass', ro*A)  
import openseespy.postprocessing.ops_vis as ops
ops.plot_model('nodes')    

########################### EIGENVALUES CALCULATION ###########################
# number of eigenvalues to calculate
eigenN = 6
# list containing lamda contaiing the first eigenN eigenvalues
lamda = op.eigen('-fullGenLapack', eigenN)    # (rad^2/s^2)
# list containing the angular frequencies of the system 
freq_Ang = []                                 # (rad/s)
for eigenvalue in lamda: 
    freq_Ang.append(eigenvalue**0.5)

op.recorder('Node', '-file', 'transDisp.out','-time', '-node',
            node_loaded, '-dof', 2, 'disp')

############################ TRANSIENT ANALYSIS ###############################
# time step
dt = 0.001                       # (s)
# final time
tEnd = 30.0                      # (s)
# number of steps
nSteps = int(tEnd / dt)
op.record()

op.setNodeDisp(8, 2, -0.07, '-commit')
op.rayleigh(0., 0., 0., 2 * 0.03 / freq_Ang[0])

op.constraints('Plain')
op.numberer('Plain')
op.system('BandGeneral')
op.test('NormDispIncr', 1e-6, nSteps)
op.algorithm('Newton')
op.integrator('Newmark', 0.5, 0.25)
op.analysis('Transient')

op.analyze(nSteps, dt)
op.wipe()

################################# PLOTTING ####################################
import numpy as np
import matplotlib.pyplot as plt
transientDisplacement = np.loadtxt("transDisp.out")
time = transientDisplacement[:, 0]
Displacement = transientDisplacement[:, 1]
plt.figure()
plt.plot(time, 1000*Displacement)
plt.xlabel('Time (s)',fontsize = 12, fontweight = 'bold')
plt.ylabel('Displacement (mm)',fontsize = 12, fontweight = 'bold')
plt.title('Free Vibration at 0.7L of the shaft',fontsize = 15, 
          fontweight = 'bold')
plt.show()

