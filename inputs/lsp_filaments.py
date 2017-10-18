######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.0035
source_loc = [0.0,-0.5,0.0087]
prop_dir = [0.0,1.0,0.0]

film_loc = [0.0,19.0,0.0087]
film_axis1 = [0.0,0.0,-0.69]
film_axis2 = [0.69,0.0,0.0]

NP = 1000000 # number of protons
ntraces = 10 # number of protons for which to track trajectories
E0 = 14.7

l_s2start = 0.499
l_prop = 0.002
nsteps = 100
spread_angle = 1.2

#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(x,y,z) that returns a tuple of all field values at any (x,y,z) point. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, ngridz, lx, ly, and lz, the number of grid cells and cell sizes you want in each dimension. Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

ngridx, ngridy, ngridz = 200, 20, 200 # REQUIRED FOR ANALYTIC GRID
grid_nthreads = 36

import numpy as np
import matplotlib.pyplot as plt

# Physical constants in cgs
q = 4.8032e-10          #statcoul
m = 1836.2*9.1094e-28   # g
c = 3.0e10              # cm/sec

# Read in LSP grid to generate physically realistic filaments
import libraries.read_xdr as lsp
FILE = lsp.flds('/g/g14/black37/LSP_files',step=369553)
(X,Y,Z,t) = FILE.get_XYZt()
dx = X[1]-X[0]
dz = Z[1]-Z[0]
lx = len(X)*dx
lz = len(Z)*dz
ly = 0.002
gridcorner = [-lx/2.0,-ly/2.0,0.0]
(By,Name,Unit) = FILE.get_VarNameUnit(name='B',xyz='y')
B_mag = max(By.flatten())*1e-4
cs = plt.contour(By,levels=[0.0])
plt.cla()
paths = cs.collections[0].get_paths()

B_fields = np.zeros((ngridx,ngridy,ngridz,3))

for i in range(ngridx):
    for j in range(ngridy):
        for k in range(ngridz):
            B_fields[i,j,k,0] += 0
            B_fields[i,j,k,1] += 0
            B_fields[i,j,k,2] += 0

# Remove any closed contours
for path in paths:
    if path.vertices[0,0] == path.vertices[-1,0] and path.vertices[0,1] == path.vertices[-1,1]:
        paths.remove(path)

# Simplify paths
for idx, path in enumerate(paths):
    paths[idx] = path.cleaned(simplify=True)

# Set up data structures for easier access to segment information, and make paths '3D' by adding random y position
for idx,path in enumerate(paths):
    new_vertices = np.zeros((len(path.vertices),3))
    new_vertices[:,0] = path.vertices[:,0]*dx-lx/2.0
    new_vertices[:,2] = path.vertices[:,1]*dz
    new_vertices[:,1] = (np.random.random()-0.5)*ly/2.0
    path.vertices = new_vertices

    path.segment_positions = (path.vertices[0:-1] + path.vertices[1:])/2.0
    path.segment_dirs = path.vertices[1:] - path.vertices[0:-1]
    path.n_segments = len(path.segment_positions)


def fields(coord):
    r = np.array(coord)
    B_vec = np.zeros(3)
    for path in paths:
        # Add contribution from each filament
        for i in range(path.n_segments):
            # Add contribution from each path segment
            dl = path.segment_dirs[i]
            rprime = r - path.segment_positions[i]
            B_vec += np.cross(dl,rprime)/np.linalg.norm(rprime)**3

    B_vec *= B_mag

    return (0.0, 0.0, 0.0, B_vec[0], B_vec[1], B_vec[2], 0.0, 0.0, 0.0)
