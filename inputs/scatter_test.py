######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.001 # Source radius in cm
source_loc = [0.0,0.0,-1.0]
prop_dir = [0.0,0.0,1.0]

film_loc = [0.0,0.0,25.0]
film_axis1 = [-1.0,0.0,0.0]
film_axis2 = [0.0,1.0,0.0]

NP = 200000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 1.0 # distance from source to z=0 (beginning of grid)
l_prop = 0.00025
nsteps = 5
spread_angle = 1 # not used

hist_maxfluence = 50.0 # optional

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

c = 3.0e10 # cm/s

ngridx, ngridy, ngridz = 1, 1, 1 # REQUIRED FOR ANALYTIC GRID

lx = 2.0 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 2.0 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = l_prop # length of simulation box in cm (REQUIRED FOR ANALYTIC GRID)

gridcorner = (-lx/2.0,-ly/2.0,0.0) # REQUIRED FOR ANALYTIC GRID

#x = np.random.normal(0,0.0001,NP)
#y = np.random.normal(0,0.0001,NP)
x = 0.0001*(np.random.random_sample(NP)-0.5)
y = 0.0001*(np.random.random_sample(NP)-0.5)
z = np.zeros(NP)
vz = np.ones(NP) * np.sqrt(2.0*E0/911.0)*c
vx = np.zeros(NP)
vy = np.zeros(NP)

fieldvals = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,79.0,5.8987546e22])

def fields(coord):
    return fieldvals


