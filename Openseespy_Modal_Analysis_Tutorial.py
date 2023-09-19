
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

# op.uniaxialMaterial('Elastic', 1)
############################# DEFINE ELEMENTS #################################
for index in range(1, Nmax):    
    op.element('elasticBeamColumn', index, *[index, index+1], A, E, Iz,
               transTag,'-mass', ro*A)
import openseespy.postprocessing.ops_vis as ops
ops.plot_model("nodes")

########################### EIGENVALUES CALCULATION #################################
  
# lamda = op.eigen('-fullGenLapack', 6)
lamda = op.eigen('-genBandArpack', 6)
# op.record()
freq_Ang = []                # (rad/s)     
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

