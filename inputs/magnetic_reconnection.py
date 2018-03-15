######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

source_fwhm = 0.0020
source_loc = [0.0,0.0,-1.0] # Location of proton source
prop_dir = [0.0,0.0,1.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0.0,0.0,27.0] # Location of center of film
film_axis1 = [-4,0.0,0.0] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,4,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 3000000 # number of protons
ntraces = 20 # number of protons for which to track trajectories
E0 = 9.5 # Initial proton energy in MeV

l_s2start = 1.0 # distance from the source to start calculating proton deflections
l_prop = 0.1 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
spread_angle = 15.0 # angle (degrees) between prop_dir and outer edge of cone of protons

hist_maxfluence = 50
plot_fluence = True # optional, default True
plot_traces = False # optional, default True
plot_quiver = False # optional, default False
save_images = False # optional, default False

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

ngridx, ngridy, ngridz = 200, 200, 100 # REQUIRED FOR ANALYTIC GRID

lx = 0.3 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 0.3 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = l_prop # REQUIRED FOR ANALYTIC GRID

gridcorner = (-lx/2.0,-ly/2.0,0.0) # REQUIRED FOR ANALYTIC GRID

grid_nthreads = 4 # Number of parallel threads used to initialize the grid (default is 1)

m = 1836.2*9.1094e-28   # g
particle_mass = 2*m # optional

Bmag = 0.4e6 # Gauss
#Emag = 1e-4*0.3333*8.0e8 # StatV/cm
Emag = 0.0
Vb = 450e5 # cm/s
t_offset = 0.3e-9
t1 = 1.3e-9
t2 = t1-t_offset

dR = 0.0100
R_outer_1 = t1*Vb
R_inner_1 = R_outer_1-dR
R_mid_1 = R_outer_1-dR/2.0

R_outer_2 = t2*Vb
R_inner_2 = R_outer_2-dR
R_mid_2 = R_outer_2-dR/2.0

x0_1,y0_1 = -0.05,-0.05
x0_2,y0_2 = 0.05,0.05

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
    Bz = 0.0
    Ex = 0.0
    Ey = 0.0
    Ez = 0.0
    N = 0
    Z = 0

    if R_inner_1 < r_1 < R_outer_1:
        theta = np.arctan2(y_1,x_1)
        R_xy = np.sqrt(x_1**2+y_1**2)
        Bx += Bmag*np.sin(theta)*(R_xy/R_outer_1)
        By += -Bmag*np.cos(theta)*(R_xy/R_outer_1)
        ex,ey,ez = Emag*np.array([x_1,y_1,z])/np.linalg.norm([x_1,y_1,z])
        Ex += ex
        Ey += ey
        Ez += ez

    if R_inner_2 < r_2 < R_outer_2:
        theta = np.arctan2(y_2,x_2)
        R_xy = np.sqrt(x_2**2+y_2**2)
        Bx += Bmag*np.sin(theta)*(R_xy/R_outer_2)
        By += -Bmag*np.cos(theta)*(R_xy/R_outer_2)
        
        ex,ey,ez = Emag*np.array([x_2,y_2,z])/np.linalg.norm([x_2,y_2,z])
        Ex += ex
        Ey += ey
        Ez += ez

    if z < 0.003:
        if (x % 0.01 < 0.003) or (y % 0.01 < 0.003):
            Z = 79.0
            N = 6.02401601e22
    return (Ex,Ey,Ez,Bx,By,Bz,0.0,Z,N)

