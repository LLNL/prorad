######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will cause AttributeError!
fformat = 'LSP' # Type of grid to load.
fname = ''
lsp_dirname = '/g/g14/black37/LSP_files'
lsp_step = 369553

lsp_ntile = 4 # Not required

r_source = 0.0 # Source radius in cm
source_loc = [0.0,-0.5,0.0087*lsp_ntile] # Location of proton source
prop_dir = [0.0,1.0,0.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0.0,19.0,0.0087*lsp_ntile] # Location of center of film
film_axis1 = [0.0,0.0,-1.38] # Vector defining direction and extent of first film axis
film_axis2 = [1.38,0.0,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 3000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 0.49 # distance from the source to start calculating proton deflections
l_prop = 0.02 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 200 # number of steps each proton takes along prop_dir
spread_angle = 6 # angle (degrees) between prop_dir and outer edge of cone of protons

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
