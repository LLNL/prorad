######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

source_fwhm = 0.0020 # Source radius in cm
source_loc = [0.0,0.0,1.35] # Location of proton source
prop_dir = [0.0,0.0,-1.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0.0,0.0,-27.0] # Location of center of film
film_axis1 = [3.0,0.0,0.0] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,3.0,0.0] # Vector defining direction and extent of second film axis (perp. to first)

# shooting along x
# source_loc = [-1.0,0.0,0.0]
# prop_dir = [1.0,0.0,0.0]
# film_loc = [27.0,0.0,0.0]
# film_axis1 = [0.0,0.0,-2.5]
# film_axis2 = [0.0,2.5,0.0]

NP = 3000000 # number of protons
ntraces = 20 # number of protons for which to track trajectories
E0 = 3.0 # Initial proton energy in MeV

l_prop = 0.1470 # distance after start (along prop_dir) through which to compute proton deflections
l_s2start = 1.275 # distance from the source to start calculating proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
spread_angle = 8.3 # angle (degrees) between prop_dir and outer edge of cone of protons

m = 1836.2*9.1094e-28   # g
particle_mass = m # optional

hist_maxfluence = 120.0

plot_traces = False

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

ngridx, ngridy, ngridz = 300, 300, 200 # REQUIRED FOR ANALYTIC GRID

lx = 0.4000 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 0.4000 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = 0.1450 # REQUIRED FOR ANALYTIC GRID

gridcorner = (-lx/2.0,-ly/2.0,-0.0050) # REQUIRED FOR ANALYTIC GRID

grid_nthreads = 4 # Number of parallel threads used to initialize the grid (default is 1)

Bmag = 0.3e6 # Gauss

Emag = 1e-4*0.3333*1.0e9 # StatV/cm

R_mid = 0.0700

def fields(coord):
    x,y,z = coord
    r = np.sqrt(x**2+y**2+z**2)
    
    Bx = 0.0
    By = 0.0
    Bz = 0.0
    Ex = 0.0
    Ey = 0.0
    Ez = 0.0
    Z = 0.0
    N = 0.0

    if z > 0:
        theta = np.arctan2(y,x)
        R_xy = np.sqrt(x**2+y**2)
        scaling = 1.0
        if R_xy < R_mid:
            scaling *= (R_xy/R_mid)
        if z > R_mid:
            scaling *= (1-(z-R_mid)/R_mid)
        mag = Bmag*2*scaling*((r-0.0250)/0.0500)**2*np.exp(-((r-0.0250)/0.0500)**2)
        Bx += mag*np.sin(theta)
        By += -mag*np.cos(theta)


    elif -0.0050 < z < 0:
        if (x % 0.0150 < 0.0050) or (y % 0.0150 < 0.0050):
            Z = 79.0
            N = 6.02401601e22
    return (Ex,Ey,Ez,Bx,By,Bz,0.0,Z,N)

