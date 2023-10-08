
@author: jarveiclos
"""
Openseespy odal analysis code designed to study a shaft with mass imbalance.
 Key parameters include a 1.5-meter shaft with a 0.5-inch diameter, a 40 kg 
 mass, a density of 7850 kg/m³, and a Young's Modulus of 2.1 GPa. Explore the 
 dynamic behavior and vibration modes of the system using this code.
"""
import openseespy.opensees as op

op.wipe()
# General model definition (2 simensions and 3 deegrees of freedom)
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

########################### EIGENVALUES CALCULATION #################################
# number of eigenvalues to calculate
eigenN = 6
# list containing lamda contaiing the first eigenN eigenvalues
lamda = op.eigen('-fullGenLapack', eigenN)    # (rad^2/s^2)
# list containing the angular frequencies of the system 
freq_Ang = []                                 # (rad/s)
for eigenvalue in lamda: 
    freq_Ang.append(eigenvalue**0.5)


# import openseespy.postprocessing.ops_vis as ops
import openseespy.postprocessing.Get_Rendering as ops

ops.plot_modeshape(1, 300)
ops.plot_modeshape(2, 300)
ops.plot_modeshape(3, 300)
ops.plot_modeshape(4, 300)
ops.plot_modeshape(5, 300)
ops.plot_modeshape(6, 300)

