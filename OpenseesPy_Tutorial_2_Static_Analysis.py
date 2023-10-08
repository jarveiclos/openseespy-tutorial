# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 06:34:07 2023

@author: carlo
"""

import openseespy.opensees as op

op.wipe()
# General model definition (2 dimensions and 3 deegrees of freedom)
op.model('Basic', '-ndm', 2, '-ndf', 3)

########################### PHYSICAL PROPERTIES ###############################
L = 1.5                  # (m)
D = 0.5 * 0.0254         # (m)
m = 40                   # (kg)
import math
pi = math.acos(-1.0)
E = 2.1e11                # (Pa)
ro = 7850                 # (kg/m^3)
A = 0.25 * pi * D**2      # (m^2)
Iz = pi * D**4 / 64       # (m^4)

############################### DEFINE NODES ##################################
Nmax = 11
for i in range(1, Nmax+1):
    op.node(i, (i - 1) * L / (Nmax - 1), 0)

############################## DOF CONSTRAINTS ################################
op.fix(   1, 1, 1, 0)
op.fix(Nmax, 1, 1, 0)

##################### ADDITION OF CONCENTRATION MASS ##########################
node_loaded = 8
op.mass(node_loaded, 0.0, m, 0.0)

transTag = 1              
op.geomTransf('Linear', transTag)

############################# DEFINE ELEMENTS #################################
for index in range(1, Nmax):    
    op.element('elasticBeamColumn', index, *[index, index+1], A, E, Iz,
               transTag,'-mass', ro*A)
import openseespy.postprocessing.ops_vis as ops
ops.plot_model("nodes")   

########################## CALCULATION OF EIGENVALUES #########################
lamda = op.eigen('-fullGenLapack', 6)
freq_Ang = []                # (rad/s)     
for eigenvalue in lamda:
    freq_Ang.append(eigenvalue**0.5)

######################### RECORDING OF EIGGENVALUES ###########################    
# # recorder eigenvalue 1
# op.recorder('Node', '-file', 'mode1VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 1')
# # recorder eigenvalue 2
# op.recorder('Node', '-file', 'mode2VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 2')
# # recorder eigenvalue 3
# op.recorder('Node', '-file', 'mode3VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 3')
# # recorder eigenvalue 4
# op.recorder('Node', '-file', 'mode4VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 4')
# # recorder eigenvalue 5
# op.recorder('Node', '-file', 'mode5VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 5')    
# # recorder eigenvalue 6
# op.recorder('Node', '-file', 'mode6VMA.out', '-nodeRange', 1 , Nmax, '-dof', 
#             2, 'eigen 6')  

# op.record()
############################# ANALYSIS DEFINITION  ############################# 
tsTag = 1
op.timeSeries('Constant', tsTag)
patternTag = 1
op.pattern('Plain', patternTag, tsTag)

op.eleLoad('-ele','-range', 1, Nmax-1, '-type', '-beamUniform',-ro*A*10)
op.load(8, 0.0, -m*10, 0.0)     # (N)

op.constraints('Plain')
op.numberer('RCM')
op.system('BandGeneral')
op.test('NormDispIncr', 1e-6, 100)
op.algorithm('Linear')
op.integrator('LoadControl', 0.1)
op.analysis('Static')

op.analyze(10)

################################# PLOTTING ####################################
import numpy as np
import matplotlib.pyplot as plt 
length = np.linspace(0, L, Nmax)

############################## SHAPE MODES ####################################
plt.figure()
for i in range(6):
    eigenMode = np.loadtxt("mode%sVMA.out"%(i+1))
    plt.plot(length, eigenMode)
plt.legend(['mode 1','mode 2','mode 3','mode 4','mode 5', 'mode 6'])
plt.show()

############################### STATIC RESPONSE ###############################
plt.figure()
ops.plot_defo()

staticDisplacement = []
for i in range(1, Nmax + 1):
    disp = op.nodeDisp(i)
    staticDisplacement.append(disp[1])
plt.figure()
plt.plot(length, staticDisplacement, marker='o')
plt.xlabel('length (m)')
plt.ylabel('y Displacement (m)')
plt.show()    

op.wipe()

