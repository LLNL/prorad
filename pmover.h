/**************************************************************************
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by M. Black amd S. C. Wilks, LLNL
 * LLNL-CODE-739358
 * All rights reserved.
 *
 * This file is part of prorad.   For details, see https://github/LLNL/prorad.
 * Please also read this link:  https://github/LLNL/prorad/AdditionalBSDNotice.
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 *
 *   *  Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the disclaimer below.
 *   *  Redistributions in binary form must reproduce the above copyright 
 *      notice, this list of conditions and the disclaimer (as noted below) 
 *      in the documentation and/or other materials provided with the 
 *      distribution.
 *   *  Neither the name of the LLNS/LLNL nor the names of its contributorsa
 *      may be used to endorse or promote products derived from this softwarea
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, 
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 ***************************************************************************/







#define EX_IDX 0
#define EY_IDX 1
#define EZ_IDX 2
#define BX_IDX 3
#define BY_IDX 4
#define BZ_IDX 5
#define NU_IDX 6
#define MATZ_IDX 7
#define DENS_IDX 8
#define NUM_FIELDS 9

#define PROT_MASS 1.6726219e-24
#define PROT_CHARGE 4.8032e-10
#define C_LIGHT 2.99792458e10
#define M_PI 3.14159265358979323846

double randNorm(void);

void scatter(double *vx, double *vy, double *vz, double Z, double N, double ds);

void pmover(int NP, int ns, double qm, double ds, double *x, double *y, double *z, 
		    double *vx, double *vy, double *vz, double *prop_dir, int ngridx, int ngridy, int ngridz,
	       	double *gridvals_flat, double *grid_offset, double *grid_spacing,
		    int ntraces, double *traces_flat, int cyl_coords);

double* nn_interp(double x, double y, double z, int ngridx, int ngridy, int ngridz,
		  double gridvals[][ngridy][ngridz][NUM_FIELDS], double *grid_offset, double *grid_spacing, int cyl_coords);
