######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

#r_source = 0.0020 # Source size in cm
source_fwhm = 0.0040
source_loc = [0,0,-1] # Location of proton source
prop_dir = [0,0,1] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0,0,27] # Location of center of film
film_axis1 = [-5,0,0] # Vector defining direction and extent of first film axis
film_axis2 = [0,5,0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 3000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 0.81 # distance from the source to start calculating proton deflections
l_prop = 0.382 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 200 # number of steps each proton takes along prop_dir
spread_angle = 15 # angle (degrees) between prop_dir and outer edge of cone of protons

hist_maxfluence = 130.0 # optional
plot_fluence = True # optional, default True
plot_traces = False # optional, default True
plot_quiver = False # optional, default False
save_images = False # optional, default False

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(x,y,z) that returns a tuple of all field values at any (x,y,z) point. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, ngridz, lx, ly, and lz, the number of grid cells and cell sizes you want in each dimension. Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

grid_nthreads = 4

cyl_coords = True

ngridx, ngridy, ngridz = 150, 450, 100 # REQUIRED FOR ANALYTIC GRID
lx,ly,lz = 0.1230, 2*np.pi, 0.38 # REQUIRED FOR ANALYTIC GRID
gridcorner = (0, 0, -0.19) # REQUIRED FOR ANALYTIC GRID

rLEH = 0.1200   # 100% LEH 2.4 mm diameter
rHOL = 0.1230   # Outer radius of the hohlraum from PRL
rCAP = 0.0250   # Capsule radius
r_AMBEFouter = 0.1100   # Outer diameter of the ambipolar field
r_AMBEFinner = 0.1050   # Inner diameter of the ambipolar field
lCAPEthickness = 0.0040  # Thickness of capsule E-field
dCAPEthickness = 0.0020  # Distance of E-field from the capsule radius

# Secondary quantities
zCAPleft = gridcorner[2] + 0.5*lz - rCAP
zCAPright = gridcorner[2] + 0.5*lz + rCAP
rCAPEinner = rCAP + dCAPEthickness
rCAPEouter = rCAPEinner + lCAPEthickness
zAMBEFleft = gridcorner[2] + 0.0*lz
zAMBEFright = gridcorner[2] + 1.0*lz

# Physical constants in cgs
q = 4.8032e-10          #statcoul
m = 1836.2*9.1094e-28   # g
c = 3.0e10              # cm/sec

# particle_charge = q # optional
particle_mass = m # optional

periodicity = 5

def fields(coord):
    r,theta,z = coord
    return (Ex(r,theta,z), Ey(r,theta,z), Ez(r,theta,z), Bx(r,theta,z), By(r,theta,z), Bz(r,theta,z), nu(r,theta,z), Z(r,theta,z), N(r,theta,z))

def Ex(r,theta,z):
    r2 = r*r
    Ex = 0.0
    if (r_AMBEFinner < r < r_AMBEFouter):
        if zAMBEFleft < z < zAMBEFright:
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ex_in_V_per_m = 7.9e8*(-np.cos(theta)*(2.0-np.cos(periodicity*theta)) \
                                   -np.sin(theta)*np.sin(periodicity*theta))
    #        - 2.0*np.sin(periodicity*theta))
            Ex = Ex_in_V_per_m*(1.0e-4)*0.3333
    if (rCAPEinner < r < rCAPEouter):
        if (zCAPleft < z < zCAPright):
            Ex_in_V_per_m = 5.6e8*np.cos(theta)  # V/m
            Ex = Ex_in_V_per_m*(1.0e-4)*0.3333   
    return Ex

def Ey(r,theta,z):
    Ey = 0.0
    if (r_AMBEFinner < r < r_AMBEFouter):
        if zAMBEFleft < z < zAMBEFright:
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ey_in_V_per_m = 7.9e8*(-np.sin(theta)*(2.0-np.cos(periodicity*theta)) \
                                   +np.cos(theta)*np.sin(periodicity*theta))
    #        - 2.0*np.sin(periodicity*theta))
            Ey = Ey_in_V_per_m*(1.0e-4)*0.3333
    if (rCAPEinner < r < rCAPEouter):
        if (zCAPleft < z < zCAPright):
            Ey_in_V_per_m = 5.6e8*np.sin(theta)  # V/m
            Ey = Ey_in_V_per_m*(1.0e-4)*0.3333
    return Ey

def Ez(r,theta,z):
    Ez = 0.0
    return Ez

def Bx(r,theta,z):
    Bx = 0.0
    return Bx

def By(r,theta,z):
    By = 0.0
    return By

def Bz(r,theta,z):
    Bz = 0.0
    return Bz

def nu(r,theta,z):
    nu = 0
    return nu

def Z(r,theta,z):
    Z = 0
    if rLEH < r < rHOL:
        Z = 79.0
    return Z

def N(r,theta,z):
    N = 0
    if rLEH < r < rHOL:
        N = 6.02401601e22
    return N
