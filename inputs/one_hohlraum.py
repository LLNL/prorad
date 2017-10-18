######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

NP = 3000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
nz_prot = 200 # number of proton z steps
l_s2grid = 0.81 # distance from source to z=0 (beginning of grid)
l_grid2film = 26.81 # distance from z=lz (end of grid) to film
film_lx = 10.0 # width of film in cm
film_ly = 10.0 # height of film in cm
Elong = 14.7 # Initial proton energy in MeV
r_source = 0.001 # Source radius in cm

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(x,y,z) that returns a tuple of all field values at any (x,y,z) point. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, ngridz, lx, ly, and lz, the number of grid cells and cell sizes you want in each dimension. Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

ngridx, ngridy, ngridz = 100, 100, 10 # REQUIRED FOR ANALYTIC GRID

lz = 0.38 # length of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lx = 0.28 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 0.28 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)   
rLEH = 0.1200   # 100% LEH 2.4 mm diameter
rHOL = 0.1230   # Outer radius of the hohlraum from PRL
rCAP = 0.0250   # Capsule radius
r_AMBEFouter = 0.1100   # Outer diamter of the ambipolar field
r_AMBEFinner = 0.1050   # Inner diameter of the ambipolar field
lCAPEthickness = 0.0040  # Thickness of capsule E-field
dCAPEthickness = 0.0020  # Distance of E-field from the capsule radius

# Secondary quantities
zCAPleft = 0.5*lz - rCAP
zCAPright = 0.5*lz + rCAP
rCAPEinner = rCAP + dCAPEthickness
rCAPEouter = rCAPEinner + lCAPEthickness
rLEH2 = rLEH*rLEH
rCAP2 = rCAP*rCAP
rHOL2 = rHOL*rHOL
r_AMBEFouter2 = r_AMBEFouter*r_AMBEFouter
r_AMBEFinner2 = r_AMBEFinner*r_AMBEFinner
rCAPEinner2 = rCAPEinner*rCAPEinner
rCAPEouter2 = rCAPEouter*rCAPEouter

# Physical constants in cgs
q = 4.8032e-10          #statcoul
m = 1836.2*9.1094e-28   # g
c = 3.0e10              # cm/sec

periodicity = 5

def fields(coord):
    x,y,z = coord
    return (Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z), nu(x,y,z),0.0,0.0)

def Bx(x,y,z):
    r2 = x**2 + y**2
    Bx = 0.0
    if ( r_AMBEFinner2 < r2 < r_AMBEFouter2 ):
        if (0.0 <= z <= lz):
            theta = np.arctan2(y, x)
            Bx_in_V_per_m = -0.0*8.0e4*np.sin(periodicity*theta)/c  # V/m
            Bx = Bx_in_V_per_m
        elif (lz <= z <= lz + 0.05):
            theta = np.arctan2(y, x)
            Bx_in_V_per_m = 0.0*8.0e4*np.sin(periodicity*theta)/c  # V/m
            Bx = Bx_in_V_per_m
    return Bx

def By(x,y,z):
    r2 = x**2 + y**2
    By = 0.0
    if ( r_AMBEFinner2 < r2 < r_AMBEFouter2 ):
        if (0.0 <= z <= lz):
            theta = np.arctan2(y, x)
            By_in_V_per_m = 0.0*8.0e4*np.cos(periodicity*theta)/c  # V/m
            By = By_in_V_per_m
        elif (lz <= z <= lz + 0.05):
            theta = np.arctan2(y, x)
            By_in_V_per_m = -0.0*8.0e4*np.cos(periodicity*theta)/c  # V/m
            By = By_in_V_per_m
    return By

def Bz(x,y,z):
    r2 = x**2 + y**2
    Bz = 0.0
    if ( 0.1**2 < r2 < 0.12**2):
        theta = np.arctan2(y, x)
        Bz_in_V_per_m = 0.0/c  # V/m
        Bz = Bz_in_V_per_m*(1.0e-4)*0.3333
    return Bz

def Ex(x,y,z):
    r2 = x**2 + y**2
    Ex = 0.0
    if ( r_AMBEFinner2 < r2 < r_AMBEFouter2 ):
        theta = np.arctan2(y, x)
        # First term is ambipolar (radial) and second term is grad Pe (theta)
        Ex_in_V_per_m = 7.9e8*(-np.cos(theta)*(2.0-np.cos(periodicity*theta)) \
                               -np.sin(theta)*np.sin(periodicity*theta))
#        - 2.0*np.sin(periodicity*theta))
        Ex = Ex_in_V_per_m*(1.0e-4)*0.3333
    if (rCAPEinner2 < r2 < rCAPEouter2):
        if (zCAPleft < z < zCAPright):
            theta = np.arctan2(y, x)
            Ex_in_V_per_m = 5.6e8*np.cos(theta)  # V/m
            Ex = Ex_in_V_per_m*(1.0e-4)*0.3333   
    return Ex

def Ey(x,y,z):
    r2 = x**2 + y**2
    Ey = 0.0
    if (r_AMBEFinner2 < r2 <r_AMBEFouter2):
        theta = np.arctan2(y, x)
        # First term is ambipolar (radial) and second term is grad Pe (theta)
        Ey_in_V_per_m = 7.9e8*(-np.sin(theta)*(2.0-np.cos(periodicity*theta)) \
                               +np.cos(theta)*np.sin(periodicity*theta))
#        - 2.0*np.sin(periodicity*theta))
        Ey = Ey_in_V_per_m*(1.0e-4)*0.3333
    if (rCAPEinner2 < r2 < rCAPEouter2):
        if (zCAPleft < z < zCAPright):
            theta = np.arctan2(y, x)
            Ey_in_V_per_m = 5.6e8*np.sin(theta)  # V/m
            Ey = Ey_in_V_per_m*(1.0e-4)*0.3333 
    return Ey

def Ez(x,y,z):
    r2 = x**2 + y**2
    Ez = 0.0
    if (0.0144 < r2 < 0.0195):
        theta = np.arctan2(y, x)
        Ez_in_V_per_m = -0.0*1.0e10   # V/m
        Ez = Ez_in_V_per_m*(1.0e-4)/3.0
    return Ez

def nu(x,y,z):
    r2 = x**2 + y**2
    nu1 = 0
    if (rLEH2 < r2 < rHOL2):
        nu1 = 0.003
    if (r2 < rCAP2):
        if (zCAPleft < z < zCAPright):
            nu1 = r2*0.003/rCAP2
    return nu1
