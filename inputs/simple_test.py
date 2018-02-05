######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

source_fwhm = 0.0
source_loc = [0,0,-1] # Location of proton source
prop_dir = [0,0,1] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0,0,5] # Location of center of film
film_axis1 = [-4,0,0] # Vector defining direction and extent of first film axis
film_axis2 = [0,4,0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 400000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 0.5 # distance from the source to start calculating proton deflections
l_prop = 1.0 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 100 # number of steps each proton takes along prop_dir
spread_angle = 3 # angle (degrees) between prop_dir and outer edge of cone of protons

m = 1836.2*9.1094e-28   # g
particle_mass = 1*m # optional
plot_traces = False

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################

ngridx, ngridy, ngridz = 1,1,1 # REQUIRED FOR ANALYTIC GRID
lx,ly,lz = 1.0,1.0,1.0 # REQUIRED FOR ANALYTIC GRID
gridcorner = (-0.5,-0.5,-0.5) # REQUIRED FOR ANALYTIC GRID

def fields(coord): # REQUIRED FOR ANALYTIC GRID
    return (0.0,0.0,0.0,1e5,1e5,0.0,0.0,0.0,0.0)
    # return (3e4,3e4,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
