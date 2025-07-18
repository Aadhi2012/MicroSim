##Geometrical dimensions of the simulation domain
DIMENSION = 3;
MESH_X = 200;
MESH_Y = 400;
MESH_Z = 200;
##Discretization, space and time
DELTA_X = 12.0;
DELTA_Y = 12.0;
DELTA_Z = 12.0;
DELTA_t = 3.0;
##Number of phases and composition
NUMPHASES = 2;
NUMCOMPONENTS = 2;
#Running and saving information
NTIMESTEPS = 10000000000;
NSMOOTH = 10;
SAVET = 10000;
STARTTIME = 0;
TRACK_PROGRESS = 200;
RESTART = 0;
numworkers = 512;
## Component and Phase names
# COMPONENTS = {Al,Cu,B};
COMPONENTS = {Al, Cu};
PHASES = {alpha, liquid};
##Material properties
##GAMMA={12, 13, 14, 23, 24...}
GAMMA = {1.0};
# Diffusivity = {Diagonal:0/1, phase, 11,22,33, 12, 13, 23...};
DIFFUSIVITY = {1, 0, 0};
DIFFUSIVITY = {1, 1, 1};
##Gas constant and molar volume
R = 1.0;
V = 1.0;
##Elasticity related parameters
#EIGEN_STRAIN = {0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
#EIGEN_STRAIN = {1, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
#VOIGT_ISOTROPIC = {0, 270, 187.5, 125.0};
#VOIGT_ISOTROPIC = {1, 270, 187.5, 125.0};
ELASTICITY = 0;
EIGEN_STRAIN = {0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0};
EIGEN_STRAIN = {1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
VOIGT_ISOTROPIC = {0, 270, 187.5, 125.0};
VOIGT_ISOTROPIC = {1, 270, 187.5, 125.0};
rho = 10;
damping_factor = 0.4;
max_iterations = 5;
#VOIGT_CUBIC = {phase, c11, c12, c44};
#VOIGT_TETRAGONAL = {phase, c11, c12, c13, c33, c44, c66};
##Boundary conditions
#0: Free, 1: Neumann, 2: Dirichlet, 3: Periodic, 4: Complex
#Boundary = {phase, X+, X-, Y+, Y-, Z+, Z-}
BOUNDARY = {phi, 3, 3, 1, 1, 3, 3};
BOUNDARY = {mu, 3, 3, 1, 1, 3, 3};
BOUNDARY = {u, 3, 3, 1, 1, 3, 3};
BOUNDARY = {c, 3, 3, 1, 1, 3, 3};
BOUNDARY = {T, 3, 3, 1, 1, 3, 3};
BOUNDARY = {F, 3, 3, 1, 1, 3, 3};
# Boundary = {phi, 1, 1, 0};
# Boundary = {"u", 3, 3, 2, 2};
#Boundary_value = {Value X+, Value X-, Value Y+, Value Y-, Value Z+, Value Z-}
BOUNDARY_VALUE = {phi, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {mu, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {c, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {T, 0, 0, 0, 0, 0, 0};
BOUNDARY_VALUE = {F, 0, 0.2, 0, 0, 0, 0};
##Type of simulation
ISOTHERMAL = 0;
BINARY = 1;
#TERNARY
DILUTE = 0;
T = 0.92;
##FILEWRITING and OUTPUTTING TO SCREEN
## WRITEFORMAT ASCII/BINARY/HDF5(Only for MPI)
##TRACK_PROGRESS: interval of writing out the progress of the simulation to stdout. 
WRITEFORMAT = BINARY;
WRITEHDF5 = 1;
##Model-specific parameters: Grand-potential model
##Phase-field parameters; epsilon:interface width; it is not the gradient energy coefficient
epsilon = 48.0;
tau = 1.31;
Tau = {0.28};
##Anisotropy functions
##Anisotropy mode, FUNCTION_ANISOTROPY=0 is isotropic
Function_anisotropy = 1;
Anisotropy_type = 4; 
dab = {0.02};
#Rotation_matrix = {0, 1, Euler_x(ang), Euler_y(ang), Euler_z(ang)};
Rotation_matrix = {0, 1, 0, 0, 0};
##Potential function
Function_W = 1;
Gamma_abc = {};
#Shifting of domain for infinite domain simulations
Shift = 1;
Shiftj = 200;
#Writing of composition fields along with the chemical potential fields
Writecomposition = 1;
#Noise
Noise_phasefield = 1;
Amp_Noise_Phase = 0.003;
##Temperature
Equilibrium_temperature = 1.0;
Filling_temperature = 1.0;
#TEMPGRADY={BASETEMP, DELTAT, DISTANCE, OFFSET, VELOCITY}
#Tempgrady = {0.96, 0.06, 800.0, 0, 0.016};
Tempgrady = {0.88, 0.12, 19200, 0, 0.010};
##Function_F
Function_F = 1;
A = {0, 1};
A = {1, 1};
ceq = {0, 0, 0.78125};
ceq = {0, 1, 0.5};
ceq = {1, 1, 0.5};
ceq = {1, 0, 0.5};
cfill = {0, 0, 0.78125};
cfill = {0, 1, 0.78125};
cfill = {1, 1, 0.78125};
cfill = {1, 0, 0.78125};
slopes = {0, 0, 0.45};
slopes = {0, 1, 0.45};
slopes = {1, 0, 0.45};
slopes = {1, 1, 0.45};
# for LBM
LBM = 1;
LBM_SAVE_FREQ = 500000000;
LBM_RESTART = 1;
NU_LBM = 0.025;
nu = 0.97;
rho_LBM = {1.2, 1.0};
beta_c  = {1.2};
gy = 0.01;
W0 = 0.8;
dt = 0.027;

