fformat = 'analytic_grid'
fname = ''
#fileout = 'hemisphere_7.45MeV.csv'

r_source = 0.0
source_loc = [-0.75, 0.0, 0.07]
prop_dir = [1.0, 0.0, 0.0]

film_loc = [12.0, 0.0, 0.07]
film_axis1 = [0.0, 0.0, 1.75]
film_axis2 = [0.0, 1.75, 0.0]

NP = 2000000
ntraces = 10
#E0 = 13.8
#E0 = 11.3
E0 = 7.45

l_s2start = 0.65
l_prop = 0.2
nsteps = 400
spread_angle = 11.5

###################


import numpy as np

ngridx, ngridy, ngridz = 200,200,100
grid_nthreads = 16

R1 = 0.0935
dR = 0.0025
#dR = 0.00016
R2 = R1+dR

R1_2 = R1**2
R2_2 = R2**2

E_SItoCGS = 1e-4*0.3333
Voltage = 5000
#Voltage = 350.0
Emag = (Voltage/(dR/100.0))*E_SItoCGS


lx = 2*R2
ly = 2*R2
lz = R2

gridcorner = (-lx/2.0, -ly/2.0, 0.0)

zeros = np.zeros(9)

def fields(coord):
    x,y,z = coord
    r_2 = (x**2+y**2+z**2)
    if R1_2 < r_2 < R2_2:
        
        # Assumes hemisphere is centered at origin!
        Edir = np.array(coord)/np.linalg.norm(coord)
        Evec = Edir*Emag
        return (Evec[0], Evec[1], Evec[2], 0.0,0.0,0.0, 0.0,0.0,0.0)
    else:
        return zeros


