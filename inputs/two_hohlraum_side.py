from two_hohlraum import *

source_loc = [-4.4,0.0,0.125]
prop_dir = [1.0,0.0,0.0] 

film_loc = [25.6,0.0,.125] # Location of center of film
film_axis1 = [0.0,0.0,4.4] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,4.4,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 4000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 3.85
l_prop = 1.1 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 400 # number of steps each proton takes along prop_dir
spread_angle = 8 # angle (degrees) between prop_dir and outer edge of cone of protons

ngridx, ngridy, ngridz = 400, 400, 50
