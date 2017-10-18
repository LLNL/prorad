######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.0035 # Source radius in cm
source_loc = [0.0,0.0,-1.0] # Location of proton source
prop_dir = [0.0,0.0,1.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0.0,0.0,27.0] # Location of center of film
film_axis1 = [4.4,0.0,0.0] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,4.4,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 5000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 15.0 # Initial proton energy in MeV

l_s2start = 1.0 # distance from the source to start calculating proton deflections
l_prop = 0.1 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
spread_angle = 10.0 # angle (degrees) between prop_dir and outer edge of cone of protons


#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

ngridx, ngridy, ngridz = 300, 300, 100 # REQUIRED FOR ANALYTIC GRID

lx = 0.3 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 0.3 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = 0.1

gridcorner = (-lx/2.0,-ly/2.0,0.0) # REQUIRED FOR ANALYTIC GRID

grid_nthreads = 16 # Number of parallel threads used to initialize the grid (default is 1)

B0 = 0.4e6 # Gauss
Vb = 450e5 # cm/s
t_offset = 0.3e-9
t1 = 1.4e-9
t2 = t1-t_offset

dR = 0.0100
R_outer_1 = t1*Vb
R_inner_1 = R_outer_1-dR

R_outer_2 = t2*Vb
R_inner_2 = R_outer_2-dR

x0_1,y0_1 = 0.05,-0.05
x0_2,y0_2 = -0.05,0.05

def fields(coord):
    x,y,z = coord
    x_1 = x-x0_1
    y_1 = y-y0_1
    x_2 = x-x0_2
    y_2 = y-y0_2
    r_1 = np.sqrt(x_1**2+y_1**2+z**2)
    r_2 = np.sqrt(x_2**2+y_2**2+z**2)
    
    Bx = 0.0
    By = 0.0

    if R_inner_1 < r_1 < R_outer_1:
        theta = np.arctan2(y_1,x_1)
        R_xy = np.sqrt(x_1**2+y_1**2)
        Bx += B0*np.sin(theta)*(R_xy/R_outer_1)
        By += -B0*np.cos(theta)*(R_xy/R_outer_1)
    if R_inner_2 < r_2 < R_outer_2:
        theta = np.arctan2(y_2,x_2)
        R_xy = np.sqrt(x_2**2+y_2**2)
        Bx += B0*np.sin(theta)*(R_xy/R_outer_2)
        By += -B0*np.cos(theta)*(R_xy/R_outer_2)

    return (0.0,0.0,0.0,Bx,By,0.0,0.0,0.0,0.0)

