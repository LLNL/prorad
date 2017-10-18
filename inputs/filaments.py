######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.0035 # Source radius in cm
source_loc = [0.0,0.0,-4.4] # Location of proton source
prop_dir = [0.0,0.0,1.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [-0.2,0.0,25.6] # Location of center of film
film_axis1 = [-1.9,0.0,0.0] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,1.9,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 3000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 4.4 # distance from the source to start calculating proton deflections
l_prop = 0.55 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
spread_angle = 6.5 # angle (degrees) between prop_dir and outer edge of cone of protons


#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

ngridx, ngridy, ngridz = 1, 300, 300 # REQUIRED FOR ANALYTIC GRID
lx = 0.55 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = 0.55 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = 0.55 # length of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
gridcorner = (-lx/2.0,-ly/2.0,0.0) # REQUIRED FOR ANALYTIC GRID

grid_nthreads = 36

import numpy as np

# Physical constants in cgs
q = 4.8032e-10          #statcoul
m = 1836.2*9.1094e-28   # g
c = 3.0e10              # cm/sec

Rot1 = np.matrix([[0,-1],[1,0]])
Rot2 = np.matrix([[0,1],[-1,0]])

nfil = 15
Rot_choices = np.random.choice(2,nfil)
Rots = []
for i in range(nfil):
    if Rot_choices[i]:
        Rots.append(Rot1)
    else:
        Rots.append(Rot2)
fil_y_coords = 0.9*(np.random.random_sample(nfil)-0.5)*ly
fil_z_coords = (0.9*np.random.random_sample(nfil)+0.05)*lz
fil_b = 0.4*lx
fil_a = 0.04*lz
B0 = -1.0e8*(fil_a/fil_b)*np.sqrt(E0)/l_s2start

def fields(coord):
    x,y,z = coord
    Bi = np.array([0.0,0.0])
    for i in range(nfil):
        Bi += B_vec(x,y,z,0.0,fil_y_coords[i],fil_z_coords[i],i)

    return (0.0, 0.0, 0.0, 0.0, Bi[0], Bi[1], 0.0, 0.0, 0.0)

def B_vec(x,y,z,x0,y0,z0,i):
    r_vec = np.array([y-y0,z-z0])
    r_mag = np.linalg.norm(r_vec)
    B_mag = B0*(r_mag/fil_a)*np.exp(-(r_mag/fil_a)**2-((x-x0)/fil_b)**2)

    r_hat = r_vec/r_mag
    return B_mag*np.array(Rots[i].dot(r_hat)).flatten()
