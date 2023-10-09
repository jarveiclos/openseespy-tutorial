
@author: jarveiclos
"""
Openseespy Modal analysis code designed to study a shaft with mass imbalance.
 Key parameters include a 1.5-meter shaft with a 0.5-inch diameter, a 40 kg 
 mass, a density of 7850 kg/m³, and a Young's Modulus of 2.1 GPa. Explore the 
 dynamic behavior and vibration modes of the system using this code.
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

########################### EIGENVALUES CALCULATION ###########################
# number of eigenvalues to calculate
eigenN = 6
# list containing lamda contaiing the first eigenN eigenvalues
lamda = op.eigen('-fullGenLapack', eigenN)    # (rad^2/s^2)
# list containing the angular frequencies of the system 
freq_Ang = []                                 # (rad/s)
for eigenvalue in lamda: 
    freq_Ang.append(eigenvalue**0.5)

########################### TRANSIENT ANALYSIS ################################    
######################## DYNAMIC LOAD DEFINITION ##############################
import numpy as np

# time step
dt = 0.001                       # (s)
# final time
tEnd = 30.0                      # (s)
# number of steps
nSteps = int(tEnd / dt)
# arrray containing all time
t = np.linspace(0, tEnd, nSteps) # (s)
# Dynamic load frequency           
omega = freq_Ang[0]              # (rad/s)
# load 1 
f1 = 0.05 * 1e-2 * m * omega**2 * np.sin(omega * t)                   # (N)
# load 2
f2 = 0.05 * 1e-2 * m * omega**2 * np.sin(0.5 * omega * t - pi / 3)    # (N)
# Resultant load 
f = f1 + f2                                                           # (N) 
np.savetxt('time.txt', t)
np.savetxt('Ampl.txt', f)
# time series Tag
tsTag = 1
op.timeSeries('Path', tsTag, '-dt', dt, '-filePath', 'Ampl.txt', '-fileTime',
                  'time.txt')
patternTag = 1
op.pattern('Plain', patternTag, tsTag)
op.load(node_loaded, 0.0, 1.0, 0.0)
# damping factor proportional to stiffness only
beta = 2 * 0.03 / freq_Ang[0]
op.rayleigh(0.0, 0.0, 0.0, beta)
   
op.recorder('Node', '-file', 'transientAn.out','-time', tsTag, '-node', 6,
                '-dof', 2, 'disp')

###########################  ANALYSIS PARAMETERS ##############################
op.constraints('Transformation')
op.numberer('Plain')
op.system('BandGeneral')
op.test('NormDispIncr', 1e-6, nSteps)
op.algorithm('Linear')
op.integrator('Newmark', 0.5, 0.25)
op.analysis('Transient')

op.record()
op.analyze(nSteps, dt)

op.wipe()
################################# PLOTTING ####################################
# list containing the recorded displacement and time
dispList = np.loadtxt('transientAn.out')
# list containing hte saved time
time = np.loadtxt('time.txt')
# list containing the saved loads
load = np.loadtxt('Ampl.txt')
############################### DYNAMIC RESPONSE ##############################
from matplotlib import pyplot as plt      
plt.plot(dispList[:,0], dispList[:,1])
plt.xlabel('time (s)')
plt.ylabel('y Displacement (m)')
################################## DYNAMIC LOAD ###############################
plt.figure()
plt.plot(time, load)
plt.xlabel('time (s)')
plt.ylabel('load (N)')
