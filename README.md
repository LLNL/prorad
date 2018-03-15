PRORAD is a general-purpose proton radiography simulation tool developed at 
NIF/LLNL. It was written by summer intern Mason Black (mason_black@hotmail.com)
and is maintained by Scott Wilks (wilks1@llnl.gov).

==================
 LLNL-CODE-739358
==================

====================
   LICENSE INFO
====================
Can be found in file https://github.com/LLNL/prorad/blob/master/LICENSE
Additional information concerning the BSD License for this software can be found 
in file: https://github.com/LLNL/prorad/blob/master/AdditionalBSDNotice

========================
    USAGE DISCLAIMER
========================

For purposes of speed, this code assumes that proton trajectories remain
approximately paraxial. What this means for practical purposes is that the 
Larmor radius of your protons, at a given energy and field strength, should 
be significantly greater than the thickness of the region over which fields of 
that strength are found. If it is not, then large errors may be introduced.

===================
    COMPILATION
===================

PRORAD is written primarily in Python 2 for ease of use, extensibility, and 
availability of libraries for reading the outputs of MHD and Rad-Hydro codes. 
However, much of the physics is handled by a compiled Python module written in 
C99, which improves performance drastically. To use this module you must first 
compile it, which is handled by a setup script. From the toplevel 'prorad' 
directory, run:

    python setup.py build

or to compile to use OpenMP (if installed on your system), run:    

    python setup.py build -omp

This will generate a file named '_pmover.so'. This can be imported by Python 
and treated as any other module. 

PRORAD can be used in two ways: interactively through a Python shell, or as a 
standalone script. In either case, the user needs to create an input file that 
defines the problem setup. The input file is itself written in Python, and must
be saved in the 'inputs' directory. See existing input files for reference.

======================
    STANDALONE USE
======================

To run as a standalone script, from the toplevel 'prorad' directory run:

    python prorad.py inputfile

where inputfile is the name of an input file within the 'inputs' directory, 
omitting the '.py'. This will load a grid, propogate protons through it, and 
plot the resulting fluence at the film plane.

=======================
    INTERACTIVE USE
=======================

To run PRORAD interactively, first start up an interactive Python session from 
the toplevel 'prorad' directory by typing 'python'.

The first step is to import prorad:

    import prorad

Then, you will need to load an input file:

    prorad.import_params('inputfile')

where 'inputfile' is whatever you named your file in the 'inputs' directory, 
omitting the '.py' (e.g. 'two_hohlraum' or 'scatter_test'). This is equivalent 
to running 'import inputs.inputfile as params' within the prorad script, in 
that all variables and functions defined in your input file will be accessible
within the 'params' namespace (e.g. 'params.l_s2grid' or 'params.fformat').

Next, you will need to generate a grid through which to propogate the protons:

    grid = prorad.load_grid()

Where the grid comes from depends on how you set 'params.fformat'. 
Current functional options are 'analytic_grid', 'FLASH', and 'HYDRA':
* 'analytic_grid' uses the user-defined function 'params.fields' from the input
  file to populate a uniformly-spaced grid. 'fields' is a function of (x,y,z)
  and returns a 9-tuple of the form (Ex,Ey,Ez,Bx,By,Bz,nu,Z,N). E is electric
  field in CGS units (StatV/cm, or V/m multiplied by 0.3333e-4), B is magnetic
  field in CGS units (Gauss), nu is stopping power (currently unused), Z is
  atomic number (used for scattering only, so set to 0.0 unless in a solid),
  and N is ion number density (cm-3, used for scattering only).
* 'FLASH' uses the `yt library <http://yt-project.org/doc/installing.html>` to 
  read in the FLASH output file specified by 'params.fname'.
* 'HYDRA' currently uses Yorick and Dave Strozzi's hyddjs tools to interpolate 
  from a 2D HYDRA dump onto a 2D x-z uniform grid, and then rotates that slice 
  around z to produce a 3D cylindrical grid. Contact David Strozzi 
  (strozzi2@llnl.gov) for access to the hyddjs tools.
* 'LSP' uses the xdrlib tool to read in LSP grids from .p4 files. Contact Drew 
  Higginson (higginson2@llnl.gov) for access.

Once you have a grid defined, you can run the simulation:

    x,y,Ep,traces = prorad.push_protons(grid)

The output tuple contains four numpy arrays. The first three are length 
params.NP and contain the final x positions, y positions, and energies of the 
protons at the film plane. The fourth is an ndarray of shape (params.ntraces, 
params.nstep+3, 3) containing the x,y,z positions of a subset of the protons at 
each point along their trajectories, including the source and film. 

Once you have run the simulation, you can analyze the final particle positions 
using whatever tools you like. To display a plot of proton fluence at the film 
plane, run:

    prorad.plot_results(grid,x,y,Ep,traces)

===================
    INPUT FILES
===================

Input files live in the 'inputs' directory. Several are included for reference.
The easiest thing to do when making a new problem setup is to copy an existing
input file and modify the values. Regardless of the type of grid being used, 
all input files must include some required parameters:

    fformat : string
        An identifier string specifying the type of grid to be used.
    fname : string
        If using a grid type that requires a data file, the absolute path to
        that file. If using an analytic grid, can leave as empty string.
    r_source : float
        Effective radius of the proton source.
    source_loc : array-like of size 3
        Location of the proton source (distances in cm).
    prop_dir : array_like of size 3
        Rather than wastefully propogate protons in all directions, the code
        will initialize protons in a cone centered about this axis. Will be
        normalized automatically.
    film_loc : array_like of size 3
        The center of the film on which final positions will be recorded.
        This, together with the two film axes, define a plane. Once the
        proton push finishes, proton trajectories are extrapolated to find
        the intersection with this plane, and the final 'x' and 'y' positions
        returned are locations on this plane relative to the two film axes.
    film_axis1 : array_like of size 3
        Vector defining the first axis of the film. Length and direction
        determine the x axis of the fluence plot. 
    film_axis2 : array_like of size 3
        Vector defining the second axis of the film. Length and direction
        determine the y axis of the fluence plot. Should be orthogonal to 
        film_axis1.
    NP : int
        Number of protons.
    ntraces : int
        Number of protons for which to record trajectories. Must be <= NP.
    E0 : float
        Initial proton energy in MeV. Determines magnitude of initial velocity.
    l_s2start : float
        Rather than being initialized at the source and propogated to the start
        of the grid, protons are instead initialized on a disc a distance
        l_s2start (in cm) along prop_dir from the source (hence the 'cone' 
        description). The distribution of positions on this disc and initial 
        velocities is set as though the protons actually did propogate here from 
        a source of the given size.
    spread_angle : float
        Half-angle (in degrees) of the cone of protons. Determines the extent
        of the disc on which protons are initialized.
    l_prop : float
        Distance (in cm) from starting plane (along prop_dir) through which to
        compute proton deflections. NOTE: this does not mean that protons only
        can move parallel to prop_dir. It only means that, during the push, 
        their net displacement projected onto prop_dir is of length l_prop. In
        addition, the displacement during each step projected onto prop_dir is 
        of equal length (l_prop/nsteps). This is forced to be the case by how
        the equations of motion are implemented in the pusher. This ensures
        that all protons start on one plane, and end on another parallel plane.
    nsteps : float
        Number of steps each proton will take during the push. The displacement
        during each step projected onto prop_dir will be an equal distance. It
        is recommended that nsteps is greater than the expected number of grid
        cells traversed along prop_dir.

There may be other required/optional parameters depending on the grid type:
    * HYDRA: 
        ngridx, ngridy, ngridz : int
            The number of radial cells, number of azimuthal (theta) cells, and 
            number of z cells. HYDRA dump is interpolated onto R-z grid and
            then rotated about z axis. 
        hyd_xrange : array-like of size 2, optional
            Min and max x values to read in from the HYDRA dump (currently the
            minimum is ignored because the cylindrical grid option requires
            the minimum radius to be 0.0).
        hyd_zrange : array-like of size 2, optional
            Min and max z values to read in from the HYDRA dump. 
        
    * analytic_grid: 
        ngridx, ngridy, ngridz : int
            The number of cells in x,y,z. 
        lx, ly, lz : float
            The spacial extend of your grid (in cm) in x,y,z.
        gridcorner : array_like of size 3
            Tuple containing the spatial position of the bottom/left/back
            corner of your grid, i.e. the position of index [0,0,0].
        fields(coord) : function, coord argument will (x,y,z) tuple.
            Returns a 9-tuple of fields at the given spatial coordinate. Will
            be called by load_grid function in prorad.py to populate a grid of
            the specified dimensions, position and extent. Try to make it fast,
            because loading large grids can take a while. The function can
            reference any other variables or functions you want to define in
            the input file. 
        grid_nthreads : int (optional)
            To speed things up, specifying a value for this will use Python's
            multiprocessing library to initialize the grid in parallel. If
            not defined the grid will be initialied in serial.
    * LSP:
        The LSP reader library currently being used requires a directory name 
        and a step number, rather than an input filename. Just leave fname as
        an empty string, and define the strings 'lsp_dirname' (full path to
        directory containing LSP dump) and 'lsp_step' (int specifying step #).
        There is also the optional parameter 'lsp_ntile' that will tile the
        grid n times in x and z (used because 2D LSP grid I tested on was small 
        but periodic).
