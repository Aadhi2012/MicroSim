#################
# Define System #
#################
```
DIMENSION = 2;
MESH_X = 90;
MESH_Y = 150;
MESH_Z = 1;
DELTA_X = 5.0;
DELTA_Y = 5.0;
DELTA_Z = 5.0;
DELTA_t = 0.3;
NUMPHASES     = 4;
NUMCOMPONENTS = 3;
```
############
# Controls #
############
```
NTIMESTEPS  = 200;
NSMOOTH     = 10;
SAVET       = 10;
RESTART     = 0;
STARTTIME   = 0;
TRACK_PROGRESS = 5;
```
#############################
# Component and Phase names #
#############################
```
COMPONENTS  = {Al, Cu, Ag};
PHASES      = {alpha, beta, gamma, liquid};
```
#################################
# Gas constant and molar volume #
################################# 
```
R = 1.0;
V = 1.0;
```
#######################
# Material properties #
#######################
# GAMMA = {12, 13, 14, 23, 24...}
```
GAMMA = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
```
# Diffusivity = {Diagonal:0/1, phase, 11,22,33, 12, 13, 23...};
```
DIFFUSIVITY = {1, 0, 0, 0, 0, 0};
DIFFUSIVITY = {1, 1, 0, 0, 0, 0};
DIFFUSIVITY = {1, 2, 0, 0, 0, 0};
DIFFUSIVITY = {1, 3, 1, 1, 0, 0};
```
#################################
# Elasticity related parameters #
#################################
ELASTICITY = 0;
```
# GRAIN_GROWTH = 0;
# rho = 10;
# damping_factor = 0.8;
# max_iterations = 5;
# EIGEN_STRAIN = {0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
# EIGEN_STRAIN = {1, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
# EIGEN_STRAIN = {2, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
# EIGEN_STRAIN = {3, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
# VOIGT_ISOTROPIC = {0, 270, 187.5, 125.0};
# VOIGT_ISOTROPIC = {1, 270, 187.5, 125.0};
# VOIGT_ISOTROPIC = {2, 270, 187.5, 125.0};
# VOIGT_ISOTROPIC = {3, 270, 187.5, 125.0};
# VOIGT_CUBIC = {phase, c11, c12, c44};
# VOIGT_TETRAGONAL = {phase, c11, c12, c13, c33, c44, c66};
```
#######################
# Boundary conditions #
#######################
```
# 0: Free, 1: Neumann, 2: Dirichlet, 3: Periodic, 4: Complex
# Boundary = {phase, X+, X-, Y+, Y-, Z+, Z-}
BOUNDARY = {phi, 3, 3, 1, 1, 3, 3};
BOUNDARY = {mu, 3, 3, 1, 1, 3, 3};
BOUNDARY = {c, 3, 3, 1, 1, 3, 3};
BOUNDARY = {T, 3, 3, 1, 1, 3, 3};
BOUNDARY = {u, 3, 3, 1, 1, 3, 3};
BOUNDARY = {F, 3, 3, 1, 1, 3, 3};
# Boundary_value = {Value X+, Value X-, Value Y+, Value Y-, Value Z+, Value Z-}
BOUNDARY_VALUE = {phi, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {mu, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {c, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {T, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {u, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {F, 0, 0, 0, 0, 0, 0};
```
#######################
## Type of simulation #
#######################
ISOTHERMAL = 0;
#DILUTE = 0;
INVERT_GSL = 0;
########################################
# FILEWRITING and OUTPUTTING TO SCREEN #
########################################
WRITEFORMAT = BINARY;
WRITEHDF5 = 1;
Writecomposition = 0;
####################################################
# Model-specific parameters: Grand-potential model #
####################################################
# epsilon : interface width
epsilon = 20.0;
tau = 1.31;
Tau = {0.837,0.837,0.837,0.837,0.837,0.837};
# Anisotropy functions 
# Anisotropy mode, FUNCTION_ANISOTROPY=0 is isotropic
Function_anisotropy = 1;
Anisotropy_type = 2; 
# dab   01,  02,  03,  12,  13,  23,
dab = {0.05, 0.0, 0.0, 0.333, 0.0, 0.0};
fab = {4.2, 0.0, 0.0, 2.0, 0.0, 0.0};
# Rotation_matrix = {0, 1, Euler_x(ang), Euler_y(ang), Euler_z(ang)};
Rotation_matrix = {0, 1, 0, 0, 0};
Rotation_matrix = {0, 2, 0, 0, 0};
Rotation_matrix = {0, 3, 0, 0, 0};
Rotation_matrix = {1, 2, 0, 0, 0};
Rotation_matrix = {1, 3, 0, 0, 0};
Rotation_matrix = {2, 3, 0, 0, 0};
# Potential function
Function_W = 1;
Gamma_abc = {7.0, 7.0, 7.0, 7.0};
######################
# Shifting of domain #
######################
Shift = 1;
Shiftj = 40;
#########
# Noise #
#########
Noise_phasefield = 0;
Amp_Noise_Phase = 0.001;
###############
# Temperature #
###############
T = 0.97;
Equilibrium_temperature = 1.0;
Filling_temperature = 1.0;
# TEMPGRADY = {BASETEMP, DELTAT, DISTANCE, OFFSET, VELOCITY}
Tempgrady = {0.92, 0.04, 800.0, 0, 0.002};
# Tempgrady = {0.92, 0.04, 800.0, 0, 0.002};
##############
# Function_F #
##############
Function_F = 1;
A = {0, 1, 1, 1, 1};
A = {1, 1, 1, 1, 1};
A = {2, 1, 1, 1, 1};
A = {3, 1, 1, 1, 1};
# liquidus 
slopes = {3, 0, 1.2, 0.0};
slopes = {0, 3, 1.2, 0.0};
slopes = {3, 1, 0.0, 1.2};
slopes = {1, 3, 0.0, 1.2};
slopes = {3, 2, -1.2, -1.2};
slopes = {2, 3, -1.2, -1.2};
# solidus
slopes = {0, 0, 1.6, 0.0};
slopes = {1, 1, 0.0, 1.6};
slopes = {2, 2, -1.6, -1.6};
ceq = {0, 0, 0.60000, 0.20000};
ceq = {0, 3, 0.33333, 0.33333};
ceq = {3, 0, 0.33333, 0.33333};
ceq = {1, 1, 0.20000, 0.60000};
ceq = {1, 3, 0.33333, 0.33333};
ceq = {3, 1, 0.33333, 0.33333};
ceq = {2, 2, 0.20000, 0.20000};
ceq = {2, 3, 0.33333, 0.33333};
ceq = {3, 2, 0.33333, 0.33333};
ceq = {3, 3, 0.33333, 0.33333};
#################
# Filling comps #
#################
cfill = {0, 0, 0.60000, 0.20000};
cfill = {0, 3, 0.33333, 0.33333};
cfill = {3, 0, 0.33333, 0.33333};
cfill = {1, 1, 0.20000, 0.60000};
cfill = {1, 3, 0.33333, 0.33333};
cfill = {3, 1, 0.33333, 0.33333};
cfill = {2, 2, 0.20000, 0.20000};
cfill = {2, 3, 0.33333, 0.33333};
cfill = {3, 2, 0.33333, 0.33333};
cfill = {3, 3, 0.33333, 0.33333};
##########
# exp bc #
##########
eidt_mode = 2;
comp_ff = {0.333, 0.333};
rad_coeff = 0.003;
step_coeff = 1e-5;
step_zero = 0;
comp_ff_rate = {0, 0};
# comp_ff_rate = {1e-12, 2e-12};
############
# LBM VARS #
############
LBM = 1;
LBM_SAVE_FREQ = 10;
LBM_RESTART = 0;
NU_LBM = 0.025;
nu = 0.97;
rho_LBM = {1.2, 1.2, 1.2, 1.0};
beta_c  = {1.0, 1.0};
gy = 0.01;
W0 = 0.8;
dt = 0.015;
# END OF FILE