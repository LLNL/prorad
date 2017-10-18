######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will cause AttributeError!
fformat = 'HYDRA' # Type of grid to load.
#fname = '/p/lscratche/kemp10/hydra/titan/v0.3/hydr02370.root' # Can leave as empty string if using fformat='analytic_grid'
fname = '/g/g14/black37/titan/v0.3/hydr02370.root' # Can leave as empty string if using fformat='analytic_grid'
#fname = '/g/g14/black37/titan/v0.3/hydr00722.root' # Can leave as empty string if using fformat='analytic_grid'
#fname = '/g/g14/black37/titanA_01/hydr00000.root' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.0 # Source radius in cm
source_loc = [-1.0,0.00,0.25] # Location of proton source
prop_dir = [1.0,0.0,0.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [1.0,-0.03,0.25] # Location of center of film
film_axis1 = [0.0,0.0,0.27] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,0.20,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 4000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 10.0 # Initial proton energy in MeV

l_s2start = 0.87 # distance from the source to start calculating proton deflections
#l_s2start = 0.01 # distance from the source to start calculating proton deflections
l_prop = 0.26 # distance after start (along prop_dir) through which to compute proton deflections
#l_prop = 1.98 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
#nsteps = 2000 # number of steps each proton takes along prop_dir
spread_angle = 9.0 # angle (degrees) between prop_dir and outer edge of cone of protons

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(p) that returns a tuple of all field values at any (x,y,z) point. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid' or 'HYDRA', you may not need to put anything here. 
"""

ngridx, ngridy, ngridz = 300, 400, 1000 # Required for analytic grid and for HYDRA

hyd_xrange = [0.0,0.12]
hyd_zrange = [-0.21,0.4]
