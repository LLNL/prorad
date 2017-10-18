######################
# REQUIRED VARIABLES #
######################
# Can change values freely, but changing the names of required variables will give import errors!
fformat = 'analytic_grid' # Type of grid to load.
fname = '' # Can leave as empty string if using fformat='analytic_grid'

r_source = 0.0035 # Source radius in cm
source_loc = [0.0,0.0,-4.4] # Location of proton source
prop_dir = [0.0,0.0,1.0] # Direction of propogation (cone of protons centered around this axis). Need not be normalized

film_loc = [0.0,0.0,25.6] # Location of center of film
film_axis1 = [-4.4,0.0,0.0] # Vector defining direction and extent of first film axis
film_axis2 = [0.0,4.4,0.0] # Vector defining direction and extent of second film axis (perp. to first)

NP = 4000000 # number of protons
ntraces = 15 # number of protons for which to track trajectories
E0 = 14.7 # Initial proton energy in MeV

l_s2start = 4.4 # distance from the source to start calculating proton deflections
l_prop = 0.25 # distance after start (along prop_dir) through which to compute proton deflections
nsteps = 500 # number of steps each proton takes along prop_dir
spread_angle = 8 # angle (degrees) between prop_dir and outer edge of cone of protons


#####################################
# USER-DEFINED VARIABLES AND FIELDS #
#####################################
"""
If you are using fformat='analytic_grid', you need to define a function fields(coord) that returns a tuple of all field values at any (x,y,z) coordinate. Field values are (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) where nu is stopping power, Z is atomic number, and N is number density. If it's not a solid at that point, set Z to zero and the scattering model won't be used. You will also need to define ngridx, ngridy, and ngridz (the number of grid cells you want in each dimension), lx, ly, and lz (the length of your grid in each dimension). If a lot of your problem has zero fields, you may want to also define field_xrange_idx, field_yrange_idx, and field_zrange_idx, which are arrays that contain the indices at which you wish to look up the fields that you define (if left undefined, all indices will be used). Define any other local variable names you want to help set up your field geometry, as long as they don't conflict with the required ones above. If you are not using fformat='analytic_grid', you probably don't need to put anything here. 
"""

import numpy as np

ngridx, ngridy, ngridz = 300, 300, 50 # REQUIRED FOR ANALYTIC GRID

lh = 0.2500 # Hohlraum length
l_field_ext_Al = 0.1500 # distance that field extends out of either end of hohlraum
l_field_ext_Au = 0.1200 # distance that field extends out of either end of hohlraum
l_field_ext = max(l_field_ext_Au,l_field_ext_Al) # how much grid extends past hohlraum
#z1 = l_field_ext # beginning of hohlraum
z1 = 0.0 # beginning of hohlraum
z2 = z1 + lh # end of hohlraum
rLEH = 0.5400/2.0   # LEH radius in cm
rHOL = rLEH+0.0050   # Outer hohlraum radius in cm
lx = rHOL*4.0 # width of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
ly = rHOL*4.0 # height of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
#lz = l_field_ext+lh+l_field_ext # length of simulation box in cm (REQUIRED FOR ANALYTIC GRID)
lz = lh

gridcorner = (-lx/2.0,-ly/2.0,0.0) # REQUIRED FOR ANALYTIC GRID

grid_nthreads = 36 # Number of parallel threads used to initialize the grid (default is 1)

# OPTIONAL: If you would like to only initialize parts of your grid (and leave the rest zero),
# you can speed up your initialization by defining which indices to initialize in each dimension.
#field_xrange_idx = range(int(ngridx*(lx/2.0-rHOL)/lx)-1,int(ngridx*(lx/2.0+rHOL)/lx)+1)
#field_yrange_idx = range(ngridy/2)
#field_zrange_idx = range(int((z1/lz)*ngridz),int((z2/lz)*ngridz)) #+ range(int((z2/lz)*ngridz), ngridz)

r_AMBEFouter_Au = rLEH-0.0100 # Outer field radius for Au in cm
r_AMBEFinner_Au = r_AMBEFouter_Au-0.0050 # Inner field radius for Au in cm
r_AMBEFouter_Al = rLEH-0.0400 # Outer field radius for Al in cm
r_AMBEFinner_Al = r_AMBEFouter_Al-0.0050 # Inner field radius for Al in cm

r2LEH = rLEH**2
r2HOL = rHOL**2
r2_AMBEFouter_Au = r_AMBEFouter_Au**2 
r2_AMBEFinner_Au = r_AMBEFinner_Au**2
r2_AMBEFouter_Al = r_AMBEFouter_Al**2
r2_AMBEFinner_Al = r_AMBEFinner_Al**2

q = 4.8032e-10          #statcoul
m = 1836.2*9.1094e-28   # g
c = 3.0e10              # cm/sec

periodicity = 16
Bfield_mag = 2e6
#Efield_mag = 6.9e8 # V/m (SI units)
Efield_mag = 4.0e8 # V/m (SI units)
E_SItoCGS = 1e-4*0.3333
gradPe_strength = 0.4
outerfield_multiplier = 0.0
Al_rot = 82*np.pi/180.0
Au_rot = Al_rot + np.pi
x0_Al, y0_Al = rHOL*np.cos(Al_rot), rHOL*np.sin(Al_rot)
x0_Au, y0_Au = rHOL*np.cos(Au_rot), rHOL*np.sin(Au_rot)

ext_field_curve_Au = l_field_ext_Au/r2HOL
ext_field_curve_Al = l_field_ext_Al/r2HOL

Ext_field_1_outer_Au = lambda r2: (z1-l_field_ext_Au) + l_field_ext_Au*r2/r2HOL
Ext_field_1_inner_Au = lambda r2: (z1-l_field_ext_Au) + l_field_ext_Au*r2/r2HOL + (r2HOL-r2LEH)*l_field_ext_Au/r2HOL
Ext_field_2_outer_Au = lambda r2: -l_field_ext_Au*r2/r2HOL + z2 + l_field_ext_Au
Ext_field_2_inner_Au = lambda r2: -l_field_ext_Au*r2/r2HOL + z2 + l_field_ext_Au - (r2HOL-r2LEH)*l_field_ext_Au/r2HOL
Ext_field_1_outer_Al = lambda r2: (z1-l_field_ext_Al) + l_field_ext_Al*r2/r2HOL
Ext_field_1_inner_Al = lambda r2: (z1-l_field_ext_Al) + l_field_ext_Al*r2/r2HOL + (r2HOL-r2LEH)*l_field_ext_Al/r2HOL
Ext_field_2_outer_Al = lambda r2: -l_field_ext_Al*r2/r2HOL + z2 + l_field_ext_Al
Ext_field_2_inner_Al = lambda r2: -l_field_ext_Al*r2/r2HOL + z2 + l_field_ext_Al - (r2HOL-r2LEH)*l_field_ext_Al/r2HOL

rot_mod = lambda x,y: (np.sin(periodicity*np.arctan2(y,x))+3)/4


def fields(coord):
    """
    Calculate fields at point coord, where coord is a tuple of (x,y,z).
    Returns tuple or array-like of field values (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N)
    """
    x,y,z = coord
    return (Ex(x,y,z), Ey(x,y,z), Ez(x,y,z), Bx(x,y,z), By(x,y,z), Bz(x,y,z), nu(x,y,z), Z(x,y,z), N(x,y,z))
    #return (Ex(x,y,z), Ey(x,y,z), 0.0, Bx(x,y,z), By(x,y,z), 0.0, 0.0, Z(x,y,z), N(x,y,z))



def Bx(x,y,z):
    Bx = 0.0
    
    r2Au = (x-x0_Au)**2 + (y-y0_Au)**2
    r2Al = (x-x0_Al)**2 + (y-y0_Al)**2

    if r2_AMBEFinner_Al < r2Al < r2_AMBEFouter_Al:
        theta = np.arctan2(y-y0_Al, x-x0_Al)
        if z < 0.0050:
            Bx = Bfield_mag*np.sin(theta)
        elif z > lh-0.0050: 
            Bx = -Bfield_mag*np.sin(theta)
    elif r2_AMBEFinner_Au < r2Au < r2_AMBEFouter_Au:
        theta = np.arctan2(y-y0_Au, x-x0_Au)
        if z < 0.0050:
            Bx = Bfield_mag*np.sin(theta)
        elif z > lh-0.0050: 
            Bx = -Bfield_mag*np.sin(theta)

    return Bx

def By(x,y,z):
    By = 0.0
    
    r2Au = (x-x0_Au)**2 + (y-y0_Au)**2
    r2Al = (x-x0_Al)**2 + (y-y0_Al)**2

    if r2_AMBEFinner_Al < r2Al < r2_AMBEFouter_Al:
        theta = np.arctan2(y-y0_Al, x-x0_Al)
        if z < 0.0050:
            By = -Bfield_mag*np.cos(theta)
        elif z > lh-0.0050: 
            By = Bfield_mag*np.cos(theta)
    elif r2_AMBEFinner_Au < r2Au < r2_AMBEFouter_Au:
        theta = np.arctan2(y-y0_Au, x-x0_Au)
        if z < 0.0050:
            By = -Bfield_mag*np.cos(theta)
        elif z > lh-0.0050: 
            By = Bfield_mag*np.cos(theta)

    return By

def Bz(x,y,z):
    return 0.0

def Ex(x,y,z):
    Ex = 0.0
    r2Au = (x-x0_Au)**2 + (y-y0_Au)**2
    r2Al = (x-x0_Al)**2 + (y-y0_Al)**2
        
    if z < z1:
        if r2Au < r2HOL and Ext_field_1_outer_Au(r2Au) < z < Ext_field_1_inner_Au(r2Au):
            Ex = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Au*(x-x0_Au) \
                    /np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_1_outer_Al(r2Al) < z < Ext_field_1_inner_Al(r2Al):
            Ex = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Al*(x-x0_Al) \
                    /np.sqrt(1+4*ext_field_curve_Al**2*r2Al)
    elif z > z2:
        if r2Au < r2HOL and Ext_field_2_inner_Au(r2Au) < z < Ext_field_2_outer_Au(r2Au):
            Ex = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Au*(x-x0_Au) \
                    /np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_2_inner_Al(r2Al) < z < Ext_field_2_outer_Al(r2Al):
            Ex = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Al*(x-x0_Al) \
                    /np.sqrt(1+4*ext_field_curve_Al**2*r2Al)
    else:
        if r2_AMBEFinner_Au < r2Au < r2_AMBEFouter_Au:
            theta = np.arctan2(y-y0_Au, x-x0_Au)
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ex_in_V_per_m = Efield_mag*(-np.cos(theta)*(2.0-np.cos(periodicity*theta)) \
                                        -gradPe_strength*np.sin(theta)*np.sin(periodicity*theta))
            Ex = Ex_in_V_per_m*E_SItoCGS
        elif r2_AMBEFinner_Al < r2Al < r2_AMBEFouter_Al:
            theta = np.arctan2(y-y0_Al, x-x0_Al)
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ex_in_V_per_m = Efield_mag*(-np.cos(theta)*(2.0-np.cos(periodicity*theta)) \
                                        -gradPe_strength*np.sin(theta)*np.sin(periodicity*theta))
            Ex = Ex_in_V_per_m*E_SItoCGS
        Ex *= np.sin((z-z1)*np.pi/lh)

    return Ex

def Ey(x,y,z):
    Ey = 0.0
    r2Au = (x-x0_Au)**2 + (y-y0_Au)**2
    r2Al = (x-x0_Al)**2 + (y-y0_Al)**2
    
    if z < z1:
        if r2Au < r2HOL and Ext_field_1_outer_Au(r2Au) < z < Ext_field_1_inner_Au(r2Au):
            Ey = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Au*(y-y0_Au) \
                    /np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_1_outer_Al(r2Al) < z < Ext_field_1_inner_Al(r2Al):
            Ey = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Al*(y-y0_Al) \
                    /np.sqrt(1+4*ext_field_curve_Al**2*r2Al)
    elif z > z2:
        if r2Au < r2HOL and Ext_field_2_inner_Au(r2Au) < z < Ext_field_2_outer_Au(r2Au):
            Ey = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Au*(y-y0_Au) \
                    /np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_2_inner_Al(r2Al) < z < Ext_field_2_outer_Al(r2Al):
            Ey = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS*2*ext_field_curve_Al*(y-y0_Al) \
                    /np.sqrt(1+4*ext_field_curve_Al**2*r2Al)
    else:
        if r2_AMBEFinner_Au < r2Au < r2_AMBEFouter_Au:
            theta = np.arctan2(y-y0_Au, x-x0_Au)
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ey_in_V_per_m = Efield_mag*(-np.sin(theta)*(2.0-np.cos(periodicity*theta)) \
                                        +gradPe_strength*np.cos(theta)*np.sin(periodicity*theta))
            Ey = Ey_in_V_per_m*E_SItoCGS
        elif r2_AMBEFinner_Al < r2Al < r2_AMBEFouter_Al:
            theta = np.arctan2(y-y0_Al, x-x0_Al)
            # First term is ambipolar (radial) and second term is grad Pe (theta)
            Ey_in_V_per_m = Efield_mag*(-np.sin(theta)*(2.0-np.cos(periodicity*theta)) \
                                        +gradPe_strength*np.cos(theta)*np.sin(periodicity*theta))
            Ey = Ey_in_V_per_m*E_SItoCGS
        Ey *= np.sin((z-z1)*np.pi/lh)

    return Ey

def Ez(x,y,z):
    Ez = 0.0
    r2Au = (x-x0_Au)**2 + (y-y0_Au)**2
    r2Al = (x-x0_Al)**2 + (y-y0_Al)**2
    
    if z < z1:
        if r2Au < r2HOL and Ext_field_1_outer_Au(r2Au) < z < Ext_field_1_inner_Au(r2Au):
            Ez = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS/np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_1_outer_Al(r2Al) < z < Ext_field_1_inner_Al(r2Al):
            Ez = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS/np.sqrt(1+4*ext_field_curve_Al**2*r2Al)
    elif z > z2:
        if r2Au < r2HOL and Ext_field_2_inner_Au(r2Au) < z < Ext_field_2_outer_Au(r2Au):
            Ez = rot_mod(x-x0_Au,y-y0_Au)*outerfield_multiplier*Efield_mag*E_SItoCGS/np.sqrt(1+4*ext_field_curve_Au**2*r2Au)
        elif r2Al < r2HOL and Ext_field_2_inner_Al(r2Al) < z < Ext_field_2_outer_Al(r2Al):
            Ez = rot_mod(x-x0_Al,y-y0_Al)*outerfield_multiplier*Efield_mag*E_SItoCGS/np.sqrt(1+4*ext_field_curve_Al**2*r2Al)

    return Ez

def nu(x,y,z):
    return 0.0

def Z(x,y,z):
    Zmat = 0.0
    if z1 < z < z2:
        r2Au = (x-x0_Au)**2+(y-y0_Au)**2
        r2Al = (x-x0_Al)**2+(y-y0_Al)**2
        if r2LEH < r2Au < r2HOL:
            """
            if abs(y-y0_Au) < 0.7*rHOL and z1+lh/8.0 < z < z2 - lh/8.0:
                Zmat = 0.0 # is a window
            else:
                Zmat = 79.0 # Au
            """
            Zmat = 79.0 # Au
        elif r2LEH < r2Al < r2HOL:
            """
            if abs(y-y0_Al) < 0.7*rHOL and z1+lh/8.0 < z < z2 - lh/8.0:
                Zmat = 0.0 # is a window
            else:
                Zmat = 13.0 # Al
            """
            Zmat = 13.0 # Al

    return Zmat

def N(x,y,z):
    Ndens = 0.0
    if z1 < z < z2:
        r2Au = (x-x0_Au)**2+(y-y0_Au)**2
        r2Al = (x-x0_Al)**2+(y-y0_Al)**2
        if r2LEH < r2Au < r2HOL:
            Ndens = 5.8987546e22 # Au
        elif r2LEH < r2Al < r2HOL:
            Ndens = 6.02401601e22 # Al
    return Ndens

