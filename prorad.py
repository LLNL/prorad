#!/usr/bin/env python
# *************************************************************************
# * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# * Produced at the Lawrence Livermore National Laboratory
# * Written by M. Black amd S. C. Wilks, LLNL
# * LLNL-CODE-739358
# * All rights reserved.
# *
# * This file is part of prorad.   For details, see https://github/LLNL/prorad.
# * Please also read this link:  https://github/LLNL/prorad/AdditionalBSDNotice.
# *
# * Redistribution and use in source and binary forms, with or without 
# * modification, are permitted provided that the following conditions are met:
# *
# *   *  Redistributions of source code must retain the above copyright notice, 
# *      this list of conditions and the disclaimer below.
# *   *  Redistributions in binary form must reproduce the above copyright 
# *      notice, this list of conditions and the disclaimer (as noted below) 
# *      in the documentation and/or other materials provided with the 
# *      distribution.
# *   *  Neither the name of the LLNS/LLNL nor the names of its contributorsa
# *      may be used to endorse or promote products derived from this softwarea
# *      without specific prior written permission.
# *
# * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, 
# * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# ***************************************************************************/

import _pmover
import numpy as np
import scipy as sp
import scipy.linalg as linalg
from scipy import interpolate
import matplotlib as mpl
import os

if not os.environ.has_key('DISPLAY'):
    mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import sys
import importlib
from multiprocessing import Pool
from timeit import default_timer as timer

sys.dont_write_bytecode = True

# Physical constants in cgs
proton_charge = 4.8032e-10  # statcoul
charge_SI = 1.6022e-19  # coulombs
proton_mass = 1836.2 * 9.1094e-28  # g
c = 2.9979e10  # cm/sec
c_SI = c * 1e-2
E_SItoCGS = 1e6 / c
MeV_to_ergs = 1.602e-6
MeV_to_Joules = 1.602e-13

# Number of field values stored at each grid point (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N)
NUM_FIELDS = 9


def main():
    # Load input parameters
    try:
        import_params(sys.argv[1])
    except (ImportError, IndexError):
        print(
            "ERROR: You must supply as a command line argument the name of a valid input parameter file contained in 'inputs/'")
        print("       (omit the path and '.py' extension)")
        return

    # Load in the grid
    grid = load_grid(fformat=params.fformat, fname=params.fname)
    if grid is None: quit()

    # Push the protons through the grid to obtain final locations on film plane
    film_x, film_y, Ep, traces = push_protons(grid)

    # Default plotting options
    plot_fluence = True
    plot_traces = True
    plot_quiver = False
    save_images = False

    # Overwrite plotting options with any that are specified in the input file
    try:
        plot_fluence = params.plot_fluence
    except AttributeError:
        pass
    try:
        plot_traces = params.plot_traces
    except AttributeError:
        pass
    try:
        plot_quiver = params.plot_quiver
    except AttributeError:
        pass
    try:
        save_images = params.save_images
    except AttributeError:
        pass

    cbar, fluence_fig, traces_fig, quiver_fig = plot_results(grid, film_x, film_y, Ep, traces,
                                                             plot_fluence=plot_fluence, plot_traces=plot_traces,
                                                             plot_quiver=plot_quiver, save_images=save_images)
    plt.show()

    # If the user has defined 'fileout', save final x and y proton positions to a file with that name in 'outputs' directory.
    try:
        save_results(film_x, film_y, params.fileout)
    except AttributeError:
        pass


def import_params(inputfile):
    """
    Load simulation parameters and put params module in global namespace.
    Normally, this will be called by main and passed whatever input name
    is supplied by the user on the command line. If prorad is imported
    in an interactive python shell, the user must call this functon to
    import the desired input parameter file.
    Parameters
    _________
    inputfile : string
        Name of an input file contained in the 'inputs' directory,
        omitting the '.py' or '.pyc' suffix.
    """
    globals()['params'] = importlib.import_module('.' + inputfile, 'inputs')


class Grid(object):
    """

    Stores field values on a grid in 3 dimensions.
    Parameters
    ----------
    gridvals : 4D ndarray with shape (nx,ny,nz,NUM_FIELDS)
        All 9 field values (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N) defined on a 3D grid
    gridspacings : tuple of floats (dx,dy,dz) or (dR,dTheta,dz)
        Spacing between grid points along each dimension
    gridcorner : tuple of floats (xoffset,yoffset,zoffset)
        Spatial coordinates of where the corner of the grid (gridvals[0,0,0,:]) is located
        (typically (-lx/2,-ly/2,0.0) for cartesian grid, (0.0,0.0,0.0) for cylindrical)
    gridextent : tuple of floats (lx, ly, lz)
        Length, width and depth of grid in Cartesian coordinates. For cylindrical grid,
        lx and ly should be at least the diameter.
    """

    def __init__(self, gridvals, gridspacings, gridcorner, gridextent, cyl_coords=False):
        self.vals = gridvals
        self.xoffset, self.yoffset, self.zoffset = gridcorner
        self.dx, self.dy, self.dz = gridspacings
        self.nx = len(gridvals[:, 0, 0, 0])
        self.ny = len(gridvals[0, :, 0, 0])
        self.nz = len(gridvals[0, 0, :, 0])
        self.lx, self.ly, self.lz = gridextent
        self.cyl_coords = cyl_coords


def load_grid(fformat=None, fname=None, ngridx=None, ngridy=None, ngridz=None):
    """
    Produce a Grid object using either analytic definitions or the output file of an MHD simulation.
    Parameters
    ----------
    fformat : string
        The type of input grid to be used. Options are 'analytic_grid', 'FLASH', 'HYDRA', or 'LSP'.
    fname : string, optional
        The name of the file to be read in, if fformat is something other than analytic_grid.
        If the file is not in the runtime directory, fname should include file path.
    ngridx : int, optional
        Desired number of grid points along x dimension. Must be either passed as an argument
        OR defined in the input file. Can be left undefined for FLASH grids.
    ngridy : int, optional
        Desired number of grid points along y dimension. Must be either passed as an argument
        OR defined in the input file. Can be left undefined for FLASH grids.
    ngridz : int, optional
        Desired number of grid points along z dimension. Must be either passed as an argument
        OR defined in the input file. Can be left undefined for FLASH grids.
    Returns
    -------
    grid : Grid object populated with field values, grid spacings, grid offset, and grid dimensions.
    """

    start_time = timer()

    # If grid dimensions not passed as arguments, try to get them from input file
    if ngridx is None:
        try:
            ngridx = params.ngridx
        except AttributeError:
            pass
    if ngridy is None:
        try:
            ngridy = params.ngridy
        except AttributeError:
            pass
    if ngridz is None:
        try:
            ngridz = params.ngridz
        except AttributeError:
            pass

    grid = None

    if fformat is None:
        try:
            fformat = params.fformat
        except AttributeError:
            print("ERROR: No file format supplied. Aborting.")
            return None

    if fname is None and (fformat != 'analytic_grid' and fformat != 'LSP'):
        try:
            fname = params.fname
        except AttributeError:
            print("ERROR: File name required for fformat='" + fformat + "'. Aborting.")

    if fformat == "analytic_grid":
        # Load user-defined analytic fields

        print("Generating " + str(ngridx) + "x" + str(ngridy) + "x" + str(ngridz) + " grid from analytic fields...")

        cyl_coords = False
        try:
            cyl_coords = params.cyl_coords
        except AttributeError:
            pass

        lx, ly, lz = params.lx, params.ly, params.lz
        xoffset, yoffset, zoffset = params.gridcorner

        X_coords = np.linspace(xoffset, xoffset + lx, ngridx)
        Y_coords = np.linspace(yoffset, yoffset + ly, ngridy)
        Z_coords = np.linspace(zoffset, zoffset + lz, ngridz)

        # Which grid indices to populate. Omitted indices will be left as zeros.
        # All indices are populated unless specified otherwise in the input file.
        # Leaving out indices that the user knows will have zero fields can speed things up.
        field_xrange_idx = range(ngridx)
        field_yrange_idx = range(ngridy)
        field_zrange_idx = range(ngridz)

        try:
            field_xrange_idx = params.field_xrange_idx
        except AttributeError:
            pass
        try:
            field_yrange_idx = params.field_yrange_idx
        except AttributeError:
            pass
        try:
            field_zrange_idx = params.field_zrange_idx
        except AttributeError:
            pass

        gridvals = np.zeros((ngridx, ngridy, ngridz, NUM_FIELDS))

        try:
            # If grid_nthreads is defined in input file, initialize the specified grid elements in parallel.
            # Using Python's multiprocessing library rather than mpi, so actual parallelism is currently limited
            # to the number of processors on a single node.
            p = Pool(params.grid_nthreads)
            print("Using " + str(params.grid_nthreads) + " threads to initialize grid")
            coords = np.meshgrid(X_coords[field_xrange_idx], Y_coords[field_yrange_idx], Z_coords[field_zrange_idx],
                                 indexing='ij')
            coords = itertools.product(X_coords[field_xrange_idx], Y_coords[field_yrange_idx],
                                       Z_coords[field_zrange_idx])
            initialized_vals = np.array(p.map(params.fields, coords))
            initialized_vals = initialized_vals.reshape(len(field_xrange_idx), len(field_yrange_idx),
                                                        len(field_zrange_idx), NUM_FIELDS)
            idx1, idx2, idx3 = np.meshgrid(field_xrange_idx, field_yrange_idx, field_zrange_idx, indexing='ij')
            gridvals[idx1, idx2, idx3] = initialized_vals
        except:
            print("'params.grid_nthreads' not specified. Initializing grid in serial.")
            for i in field_xrange_idx:
                for j in field_yrange_idx:
                    for k in field_zrange_idx:
                        x = X_coords[i]
                        y = Y_coords[j]
                        z = Z_coords[k]
                        gridvals[i, j, k, :] = params.fields((x, y, z))

        gridspacings = (lx / ngridx, ly / ngridy, lz / ngridz)

        grid = Grid(gridvals, gridspacings, (xoffset, yoffset, zoffset), (lx, ly, lz), cyl_coords=cyl_coords)

    elif fformat == "FLASH":
        print("Loading FLASH grid...")
        try:
            import yt
        except ImportError:
            print("ERROR: You need the yt module installed to load FLASH grids.")
            print("See instructions at http://yt-project.org/doc/installing.html")
            return None
        # Load the dataset using yt
        ds = yt.load(fname)

        # Sample the data onto a uniform grid, taking the coarsest resolution (i.e. averaging out any AMR)
        uniform_data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
        magx = uniform_data[u'magx'].in_cgs().to_ndarray()
        magy = uniform_data[u'magy'].in_cgs().to_ndarray()
        magz = uniform_data[u'magz'].in_cgs().to_ndarray()

        right_edge = ds.domain_right_edge.in_cgs().to_ndarray()
        left_edge = ds.domain_left_edge.in_cgs().to_ndarray()
        gridspacings = (right_edge - left_edge) / magx.shape

        ngridx, ngridy, ngridz = magx.shape

        gridvals = np.zeros((ngridx, ngridy, ngridz, NUM_FIELDS))

        for i in range(ngridx):
            for j in range(ngridy):
                for k in range(ngridz):
                    # TODO: Calculate electric fields too
                    gridvals[i, j, k, 3:6] = [magx[i, j, k], magy[i, j, k], magz[i, j, k]]

        grid = Grid(gridvals, gridspacings, left_edge, right_edge - left_edge)

    elif fformat == "HYDRA":
        print("Loading HYDRA grid...")
        try:
            import std_yorick as stdY
            import libraries.pydjs.hyddjs as DJH

        except ImportError:
            print("ERROR: You need to install and configure Yorick and the hyddjs tools to load HYDRA grids.")
            print("Contact Dave Strozzi (strozzi2@llnl.gov) for info.")
            return None

        H = stdY.h2open(fname)
        varL = ['x', 'y', 'z', 'Bx', 'By', 'Bz', 'p', 'eden', 'zb', 'ireg', 'tmat']
        stg = DJH.h2stgchk(H, varL)
        x = stg['x'][0]
        y = stg['y'][0]
        z = stg['z'][0]
        # zb = stg['zb'][0]

        xmin, xmax = 0.0, np.amax(x)
        zmin, zmax = np.amin(z), np.amax(z)

        try:
            xmin, xmax = params.hyd_xrange
        except AttributeError:
            pass
        try:
            zmin, zmax = params.hyd_zrange
        except AttributeError:
            pass

        R = np.linspace(0, xmax, ngridx)
        Theta = np.linspace(0, 2 * np.pi - 2 * np.pi / ngridy, ngridy)
        Z = np.linspace(zmin, zmax, ngridz)

        dz = Z[1] - Z[0]
        dR = R[1] - R[0]
        dTheta = Theta[1] - Theta[0]

        xz_vals = np.zeros((ngridx, ngridz, NUM_FIELDS))
        X_2D, Z_2D = np.meshgrid(R, Z, indexing='ij')

        # Calculate electric field from electron temperature and density gradient
        edens = DJH.h2interp(H, stg, 'eden', X_2D, Z_2D) * 1e6  # m^-3
        epres = DJH.h2interp(H, stg, 'p', X_2D, Z_2D) * 1e11  # Pa
        # etemp = DJH.h2interp(H,stg,'tmat',X_2D,Z_2D)*1e-3 # eV
        np.seterr(divide='ignore', invalid='ignore')  # Temporarily ignore errors because we could have divide by zeros
        # Exz = (etemp/proton_charge)*np.nan_to_num(np.divide(np.gradient(epres,dR,dz),epres))
        Exz = -np.nan_to_num(np.divide(np.gradient(epres, dR / 100.0, dz / 100.0), charge_SI * edens)) * E_SItoCGS
        np.seterr(divide='warn', invalid='warn')  # Restore to normal error behavior

        xz_vals[:, :, 0], xz_vals[:, :, 2] = Exz
        xz_vals[:, :, 3] = DJH.h2interp(H, stg, 'Bx', X_2D, Z_2D)
        xz_vals[:, :, 4] = DJH.h2interp(H, stg, 'By', X_2D, Z_2D)
        xz_vals[:, :, 5] = DJH.h2interp(H, stg, 'Bz', X_2D, Z_2D)
        # xz_vals[:,:,7] = DJH.h2interp(H,stg,'zb',X_2D,Z_2D)

        gridvals = np.zeros((ngridx, ngridy, ngridz, NUM_FIELDS))

        for j in range(ngridy):
            theta = Theta[j]
            gridvals[:, j, :, 0] = xz_vals[:, :, 0] * np.cos(theta)
            gridvals[:, j, :, 1] = xz_vals[:, :, 0] * np.sin(theta)
            gridvals[:, j, :, 2] = xz_vals[:, :, 2]
            gridvals[:, j, :, 3] = xz_vals[:, :, 3] * np.cos(theta) - xz_vals[:, :, 4] * np.sin(theta)
            gridvals[:, j, :, 4] = xz_vals[:, :, 3] * np.sin(theta) + xz_vals[:, :, 4] * np.cos(theta)
            gridvals[:, j, :, 5] = xz_vals[:, :, 5]

        grid = Grid(gridvals, (dR, dTheta, dz), (0.0, 0.0, 0.0), (2 * xmax, 2 * xmax, zmax - zmin), cyl_coords=True)

    elif fformat == "LSP":
        print("Loading LSP grid...")
        try:
            import libraries.read_xdr as lsp
        except ImportError:
            print("ERROR: You need the xdrlib tool to read LSP grids. Contact Drew Higginson for access.")
            return None

        FILE = lsp.flds(params.lsp_dirname, step=int(params.lsp_step))
        (X, Y, Z, t) = FILE.get_XYZt()

        (Ex, Name, Unit) = FILE.get_VarNameUnit(name='E', xyz='x')
        (Ey, Name, Unit) = FILE.get_VarNameUnit(name='E', xyz='y')
        (Ez, Name, Unit) = FILE.get_VarNameUnit(name='E', xyz='z')
        (Bx, Name, Unit) = FILE.get_VarNameUnit(name='B', xyz='x')
        (By, Name, Unit) = FILE.get_VarNameUnit(name='B', xyz='y')
        (Bz, Name, Unit) = FILE.get_VarNameUnit(name='B', xyz='z')

        ngridx = len(X)
        ngridy = len(Y)
        ngridz = len(Z)

        gridvals = np.zeros((ngridx, ngridy, ngridz, NUM_FIELDS))

        dx = X[1] - X[0]
        dz = Z[1] - Z[0]

        if ngridy == 1:
            # 2D grid

            # Extrude distance dy based on the wavelength of periodic features in x
            # (probably will not want to do this in general)
            By_fft = np.fft.fft(By[:, 0])
            x_periods = 1.0 / np.fft.fftfreq(len(By_fft), dx)
            dy = abs(x_periods[np.argmax(By_fft)])

            # gridvals[:,0,:,0] = Ex*1e5*E_SItoCGS
            # gridvals[:,0,:,1] = Ey*1e5*E_SItoCGS
            # gridvals[:,0,:,2] = Ez*1e5*E_SItoCGS
            gridvals[:, 0, :, 3] = Bx
            gridvals[:, 0, :, 4] = By
            gridvals[:, 0, :, 5] = Bz

            # If params.lsp_ntile is defined, tile the grid that number of times in x and z
            try:
                gridvals = np.tile(gridvals, (params.lsp_ntile, 1, params.lsp_ntile, 1))
            except AttributeError:
                pass

        else:
            # 3D grid
            # TODO: Test 3D LSP grid

            dy = Y[1] - Y[0]

            # gridvals[:,:,:,0] = Ex*1e5*E_SItoCGS
            # gridvals[:,:,:,1] = Ey*1e5*E_SItoCGS
            # gridvals[:,:,:,2] = Ez*1e5*E_SItoCGS
            gridvals[:, :, :, 3] = Bx
            gridvals[:, :, :, 4] = By
            gridvals[:, :, :, 5] = Bz

        lx = dx * len(gridvals[:, 0, 0, 0])
        ly = dy * len(gridvals[0, :, 0, 0])
        lz = dz * len(gridvals[0, 0, :, 0])

        grid = Grid(gridvals, (dx, dy, dz), (-lx / 2.0, -ly / 2.0, 0.0), (lx, ly, lz))

    else:
        print('"' + fformat + '"' + 'is not a recognized file format. Aborting.')
        return None

    end_time = timer()
    print("Time elapsed during grid generation: " + str(end_time - start_time) + " s")

    if grid.cyl_coords:
        print("Grid dR, dTheta, dz: " + str(grid.dx) + " cm, " + str(grid.dy) + " rad, " + str(grid.dz) + " cm")
        print("Grid nR, nTheta, nz: " + str(grid.nx) + ", " + str(grid.ny) + ", " + str(grid.nz))
        print("Grid lR, lTheta, lz: " + str(grid.lx) + " cm, " + str(grid.ly) + " rad, " + str(grid.lz) + " cm")
    else:
        print("Grid dx, dy, dz: " + str(grid.dx) + " cm, " + str(grid.dy) + " cm, " + str(grid.dz) + " cm")
        print("Grid nx, ny, nz: " + str(grid.nx) + ", " + str(grid.ny) + ", " + str(grid.nz))
        print("Grid lx, ly, lz: " + str(grid.lx) + " cm, " + str(grid.ly) + " cm, " + str(grid.lz) + " cm")

    return grid


def push_protons(grid):
    """
    Initialize arrays of proton positions, velocities, and energies, then simulate proton propogation through the grid.
    Parameters
    ----------
    grid : Grid
        Grid object produced by load_grid()
    Returns
    -------
    film_x : array of size NP
        'x' positions of protons on film plane relative to film axes (in cm)
    film_y : array of size NP
        'y' positions of protons on film plane relative to film axes (in cm)
    Ep : array of size NP
        Energies of protons at image plane (in MeV).
    traces : ndarray of shape (ntraces, nz_prot+2, 3)
        Record of particle traces for random subset of protons. First dimension is the proton being tracked,
        second dimension is the index of the step, third dimension is the coordinates (x,y,z) at each location
    """

    NP = params.NP
    ntraces = params.ntraces
    nsteps = params.nsteps

    # Array to store trajectories of a subset of the protons.
    # Three extra indices for when the protons are at the source, at the beginning of the grid,
    # and at the film. The 3 is for x,y,z coords of each proton.
    traces = np.zeros((ntraces, nsteps + 3, 3))

    mass = proton_mass
    charge = proton_charge

    try:
        mass = params.particle_mass
        print('Using user-defined particle mass.')
    except AttributeError:
        pass

    try:
        charge = params.particle_charge
        print('Using user-defined particle charge.')
    except AttributeError:
        pass

    try:
        # The user can choose to supply their own initial x,y,z positions and velocities.
        # (can supply either all of these or none of them, but not just some)
        # If these are provided, don't need to define params.source_loc, params.prop_dir,
        # params.l_s2start, params.E0, params.spread_angle, or params.r_source. Note that
        # protons should all be moving in a roughly similar direction, otherwise
        # large errors will be introduced because of how the equations of motion are
        # implemented in the pusher.
        x, y, z = params.x, params.y, params.z
        vx, vy, vz = params.vx, params.vy, params.vz

        prop_dir = np.array([np.average(vx), np.average(vy), np.average(vz)])
        prop_dir /= linalg.norm(prop_dir)

        print('Using user-defined initial positions and velocities.')
    except AttributeError:
        # If user didn't supply initial positions and velocities, assume isotropic source
        print('Initializing positions and velocities assuming isotropic proton source.')
        start_time = timer()

        # Protons "start" at source (meaning velocities are set as if they did)
        traces[:, 0, 0:3] = params.source_loc

        v0 = np.sqrt(2.0 * params.E0 * MeV_to_ergs / mass)

        # 3D position of proton source
        source_loc = np.array(params.source_loc)

        # unit vector defining propogation direction
        prop_dir = np.array(params.prop_dir) / linalg.norm(params.prop_dir)

        # Distance from source to where protons actually start their propogation
        l_s2start = params.l_s2start

        # Angle between prop_dir and the outer edge of the cone defining the initial proton spread
        phi_max = params.spread_angle * np.pi / 180.0

        # Together, theta and phi give an isotropic distribution of protons on the surface of a
        # unit sphere sector of angle phi_max. phis are polar angles from the distribution
        # center, and thetas are azimuthal angles.
        phis = np.arccos((1 - np.cos(phi_max)) * np.random.random_sample(NP) + np.cos(phi_max))
        thetas = 2 * np.pi * np.random.random_sample(NP)

        # Initialize at a starting plane as if source is at [0,0,0] pointing along z axis
        positions = l_s2start * (np.array([0, 0, 1])[:, np.newaxis] + np.tan(phis) * np.array(
            [np.cos(thetas), np.sin(thetas), np.zeros(NP)]))

        # Then rotate to correct direction and add in actual source position
        rot_axis = np.cross([0, 0, 1], prop_dir)
        if linalg.norm(rot_axis) == 0:
            rot_axis = [0, 1, 0]
        rot_axis /= linalg.norm(rot_axis)
        rot_angle = np.arccos(np.dot([0, 0, 1], prop_dir))
        # New favorite way to define a rotation matrix:
        rot_matrix = linalg.expm(np.cross(np.eye(3), rot_axis * rot_angle))
        positions = np.dot(rot_matrix, positions)
        positions += source_loc[:, np.newaxis]

        # Calculate velocity vector based on starting position
        s2pos = positions - source_loc[:, np.newaxis]
        velocities = v0 * s2pos / np.linalg.norm(s2pos, axis=0)

        x, y, z = positions
        vx, vy, vz = velocities

        try:
            # Add random noise to account for finite source.
            # NOTE: this means that protons might not necessarily align with
            # the starting end ending planes on either side of the grid, so
            # be careful if you have, for example, a very thin sheet of metal
            # at one end of your grid, that your l_prop will still push them
            # all through it.
            r_source = params.r_source
            x += (2 * np.random.random_sample(NP) - 1) * r_source
            y += (2 * np.random.random_sample(NP) - 1) * r_source
            z += (2 * np.random.random_sample(NP) - 1) * r_source
        except AttributeError:
            pass

        try:
            # Add random noise to account for finite source.
            # NOTE: this means that protons might not necessarily align with
            # the starting end ending planes on either side of the grid, so
            # be careful if you have, for example, a very thin sheet of metal
            # at one end of your grid, that your l_prop will still push them
            # all through it.
            source_fwhm = params.source_fwhm
            source_variance = (source_fwhm / 2.355) ** 2
            gaussian_noise = np.random.multivariate_normal([0, 0, 0], np.eye(3) * source_variance, NP)
            x += gaussian_noise[:, 0]
            y += gaussian_noise[:, 1]
            z += gaussian_noise[:, 2]
        except AttributeError:
            pass

        end_time = timer()
        print("Time elapsed during proton initialization: " + str(end_time - start_time) + " s")

    # Proton step size
    ds_prot = params.l_prop / nsteps

    print("Pushing " + str(NP) + " protons...")

    # Push protons using compiled module written in C
    start_time = timer()
    _pmover.pmover(mass, charge, ds_prot, grid.dx, grid.dy, grid.dz, grid.xoffset, grid.yoffset, grid.zoffset, NP,
                   nsteps, \
                   grid.nx, grid.ny, grid.nz, ntraces, grid.cyl_coords, x, y, z, vx, vy, vz, prop_dir, grid.vals,
                   traces)

    traces[:, nsteps + 1, 0] = x[0:ntraces]
    traces[:, nsteps + 1, 1] = y[0:ntraces]
    traces[:, nsteps + 1, 2] = z[0:ntraces]

    end_time = timer()
    print("Time elapsed during proton push through grid: " + str(end_time - start_time) + " s")

    # final propagation to film plane

    start_time = timer()

    film_axis1 = np.array(params.film_axis1) / linalg.norm(params.film_axis1)
    film_axis2 = np.array(params.film_axis2) / linalg.norm(params.film_axis2)

    film_perp = np.cross(film_axis1, film_axis2)  # Film normal vector
    film_loc = np.array(params.film_loc)

    positions = np.array([x, y, z])
    velocities = np.array([vx, vy, vz])
    positions += velocities * (
                np.dot((film_loc[:, np.newaxis] - positions).T, film_perp) / np.dot(velocities.T, film_perp))

    # Project final positions onto basis vectors that span the film
    film_x = np.dot((positions - film_loc[:, np.newaxis]).T, film_axis1)
    film_y = np.dot((positions - film_loc[:, np.newaxis]).T, film_axis2)

    x, y, z = positions

    end_time = timer()
    print("Time elapsed during final propogation and projection onto film plane: " + str(end_time - start_time) + " s")

    # Report final energy in MeV
    Ep = 0.5 * mass * (vx ** 2 + vy ** 2 + vz ** 2) / MeV_to_ergs

    traces[:, nsteps + 2, 0] = x[0:ntraces]
    traces[:, nsteps + 2, 1] = y[0:ntraces]
    traces[:, nsteps + 2, 2] = z[0:ntraces]

    film_x = np.nan_to_num(film_x)
    film_y = np.nan_to_num(film_y)
    Ep = np.nan_to_num(Ep)
    traces = np.nan_to_num(traces)

    return film_x, film_y, Ep, traces


def plot_results(grid, film_x, film_y, Ep, traces, plot_fluence=True, plot_quiver=False, plot_traces=True,
                 save_images=False):
    """
    Plot proton fluence, vector field slice, and/or 3D particle traces
    Parameters
    ----------
    grid : Grid
        Grid object produced by load_grid()
    film_x : array of size NP
        'x' positions of protons on film plane (returned by push_protons)
    film_y : array of size NP
        'y' positions of protons on film plane (returned by push_protons)
    Ep : array of size NP
        Energies of protons at film plane (in MeV).
    traces : ndarray of shape (ntraces, nz_prot+2, 3)
        Record of particle traces for random subset of protons. First dimension is the proton being tracked,
        second dimension is the index of the step, third dimension is the coordinates (x,y,z) at each location
    plot_fluence : boolean, optional
        Whether to plot proton fluence as a 2D histogram
    plot_quiver : boolean, optional
        Whether to plot a 2D x-y slice of the magnetic or electric fields using a quiver plot.
        Also plots material Z using imshow.
    plot_traces : boolean, optional
        Whether to plot particle traces in 3D using mpl_toolkits.mplot3d (if installed)
    save_images : boolean, optional
        Whether to save fluence and/or quiver plots in the 'outputs' directory. Will be saved as png images,
        and named according to the name of the input parameter file.
    """

    quiver_fig = None
    traces_fig = None
    fluence_fig = None
    cbar = None

    if plot_quiver:
        # TODO: Clean up this section

        # Plot material Z and field cross section
        quiver_fig = plt.figure(1)
        plt.axes().set_aspect('equal')

        # Material Z
        # if grid.cyl_coords is False:
        #     plt.imshow(grid.vals[:,:,grid.nz/2,7].T, origin='lower', cmap='Reds', interpolation='nearest') # plot material Z

        # x-y electric fields
        # quiverU = grid.vals[:,:,grid.nz/2,0]
        # quiverV = grid.vals[:,:,grid.nz/2,1]
        # x-y magnetic fields
        quiverU = grid.vals[:, :, grid.nz / 3, 3]
        quiverV = grid.vals[:, :, grid.nz / 3, 4]
        quiverX = range(len(quiverU[:, 0]))
        quiverY = range(len(quiverU[0, :]))
        quiverX, quiverY = np.meshgrid(quiverX, quiverY, indexing='ij')

        if grid.cyl_coords:
            # X-Y view
            quiverR = np.linspace(grid.dx / 2, grid.dx * grid.nx + grid.dx / 2, grid.nx)
            quiverTheta = np.linspace(0.0, 2 * np.pi - grid.dy, grid.ny)
            quiverR, quiverTheta = np.meshgrid(quiverR, quiverTheta, indexing='ij')
            quiverX = quiverR * np.cos(quiverTheta)
            quiverY = quiverR * np.sin(quiverTheta)

            # Z-X view (uncomment to activate)
            # quiverX = np.linspace(grid.zoffset, grid.zoffset+grid.lz, grid.nz)
            # quiverY = np.linspace(grid.dx/2,grid.dx*grid.nx+grid.dx/2,grid.nx)
            # quiverX, quiverY = np.meshgrid(quiverX,quiverY, indexing='ij')
            # quiverU = grid.vals[:,0,:,2]
            # quiverV = grid.vals[:,0,:,0]

        # Mask to get rid of the zero vectors
        quiverMask = ((quiverU != 0.0) | (quiverV != 0.0))
        plt.quiver(quiverX[quiverMask], quiverY[quiverMask], quiverU[quiverMask], quiverV[quiverMask])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim([min(quiverX.flatten()), max(quiverX.flatten())])
        plt.ylim([min(quiverY.flatten()), max(quiverY.flatten())])

        if save_images:
            plt.savefig('outputs/' + sys.argv[1] + '_quiver.png', bbox_inches='tight', transparent=True, dpi=200)

    if plot_traces:
        # Plot particle traces in 3D
        try:
            from mpl_toolkits.mplot3d import Axes3D
            traces_fig = plt.figure(2, figsize=(7, 7))
            ax = traces_fig.gca(projection='3d')
            for i in range(params.ntraces):
                # ax.plot(traces[i,0:params.nsteps+2,0],traces[i,0:params.nsteps+2,1],traces[i,0:params.nsteps+2,2])
                ax.plot(traces[i, :, 0], traces[i, :, 1], traces[i, :, 2])

            # Plot grid bounds.
            # NOTE: DO NOT MISTAKE GRID BOUNDS FOR WHERE YOUR HOHLRAUM IS, ESPECIALLY IN CYLINDRICAL CASE
            if grid.cyl_coords:
                # Plot 3d transparent cylinder
                import mpl_toolkits.mplot3d.art3d as art3d
                from matplotlib.patches import Circle
                radius = grid.lx

                # Cylinder top and bottom
                # cyl_bottom = Circle((0, 0), radius, color='k', alpha=0.2)
                # cyl_top = Circle((0, 0), radius, color='k', alpha=0.2)
                # ax.add_patch(cyl_bottom)
                # ax.add_patch(cyl_top)
                # art3d.pathpatch_2d_to_3d(cyl_top, z=grid.zoffset+grid.lz, zdir="z")
                # art3d.pathpatch_2d_to_3d(cyl_bottom, z=grid.zoffset, zdir="z")

                # Cylinder sides
                X, Z = np.meshgrid(np.linspace(-radius, radius, 20),
                                   np.linspace(grid.zoffset, grid.zoffset + grid.lz, 20))
                Y = np.sqrt(radius ** 2 - X ** 2)  # Pythagorean theorem
                ax.plot_surface(X, Y, Z, linewidth=1, color='k', alpha=0.2)
                ax.plot_surface(X, -Y, Z, linewidth=1, color='k', alpha=0.2)
            else:
                # Plot edges of 3d box as dotted black lines
                for s, e in itertools.combinations(np.array(list(
                        itertools.product([params.gridcorner[0], params.gridcorner[0] + grid.lx],
                                          [params.gridcorner[1], params.gridcorner[1] + grid.ly],
                                          [params.gridcorner[2], params.gridcorner[2] + grid.lz]))), 2):
                    if np.sum(np.abs(s - e)) in (grid.lx, grid.ly, grid.lz):
                        ax.plot3D(*zip(s, e), color='k', linestyle='--')

            if grid.cyl_coords:
                ax.set_xlim([grid.xoffset - grid.lx, grid.xoffset + grid.lx])
                ax.set_ylim([grid.xoffset - grid.lx, grid.xoffset + grid.lx])
                ax.set_zlim([grid.zoffset, grid.zoffset + grid.lz])
            else:
                ax.set_xlim([grid.xoffset, grid.xoffset + grid.lx])
                ax.set_ylim([grid.yoffset, grid.yoffset + grid.ly])
                ax.set_zlim([grid.zoffset, grid.zoffset + grid.lz])

            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        except ImportError:
            print("Unable to plot particle traces--module mpl_toolkits.mplot3d is not installed.")

    if plot_fluence:
        fluence_fig = plt.figure(3)
        plt.clf()
        plt.axes().set_aspect('equal')
        xmax = linalg.norm(params.film_axis1)
        ymax = linalg.norm(params.film_axis2)
        maxfluence = 150.0
        try:
            maxfluence = params.hist_maxfluence
        except AttributeError:
            pass
        myhist = plt.hist2d(film_x * 10, film_y * 10, bins=300, cmap='gray_r',
                            range=[[-xmax * 10, xmax * 10], [-ymax * 10, ymax * 10]], vmin=0.0, vmax=maxfluence)
        plt.xlabel('mm')
        plt.ylabel('mm')
        # plt.title('9.5 MeV deuteron')

        cbar = plt.colorbar(format='%05.2f')
        try:
            # Include interactive draggable colorbar, if available
            from libraries import draggable_cbar
            cbar = draggable_cbar.DraggableColorbar(cbar, myhist[3])
            cbar.connect()
        except ImportError:
            pass

        if save_images:
            plt.savefig('outputs/' + sys.argv[1] + '_fluence.png', bbox_inches='tight')

    return cbar, fluence_fig, traces_fig, quiver_fig


def save_results(x_arr, y_arr, fileout):
    """
    Save final proton positions to a csv output file
    Parameters
    ----------
    x_arr : array of size NP
        'x' positions of protons on film plane (returned by push_protons)
    y_arr : array of size NP
        'y' positions of protons on film plane (returned by push_protons)
    fileout : string
        What to name file in the output directory
    """

    if fileout == '': return
    f = open('outputs/' + fileout, 'w')
    for i in range(params.NP):
        x, y = x_arr[i], y_arr[i]
        f.write(str(x) + ', ' + str(y) + '\n')
    f.close()


if __name__ == '__main__':
    main()