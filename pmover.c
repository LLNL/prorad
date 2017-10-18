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

#include "pmover.h"
#ifdef USE_OMP
    #include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

static double zeros[NUM_FIELDS] = {0.0};

double* nn_interp(double x, double y, double z, int ngridx, int ngridy, int ngridz,
		  double gridvals[][ngridy][ngridz][NUM_FIELDS], double *grid_offset, double *grid_spacing, int cyl_coords) {
    
    if (cyl_coords) {
        // Cylindrical grid, need to interpolate accordingly
        int r_idx = round(sqrt(pow(x-grid_offset[0],2)+pow(y-grid_offset[1],2))/grid_spacing[0]);
        
        double theta = atan2(y-grid_offset[1],x-grid_offset[0]);
        theta = (theta >= 0.0) ? theta : theta+2*M_PI; // atan2 gives between -pi and pi, need between 0 and 2pi
        int theta_idx = round(theta/grid_spacing[1]);
        if (theta_idx == ngridy) { 
            theta_idx = 0; 
        }
        
        int z_idx = round((z - grid_offset[2])/grid_spacing[2]);
        
        if (r_idx < 0 || r_idx >= ngridx || theta_idx < 0 || theta_idx >= ngridy || z_idx < 0 || z_idx >= ngridz) {
            return zeros;
        } else {
            return &gridvals[r_idx][theta_idx][z_idx][0];
        }

    } else {
        // Cartesian grid
        int x_idx = (int)((x - grid_offset[0])/grid_spacing[0]);
        int y_idx = (int)((y - grid_offset[1])/grid_spacing[1]);
        int z_idx = (int)((z - grid_offset[2])/grid_spacing[2]);
        
        if (x_idx < 0 || x_idx >= ngridx || y_idx < 0 || y_idx >= ngridy || z_idx < 0 || z_idx >= ngridz) {
            return zeros;
        } else {
            return &gridvals[x_idx][y_idx][z_idx][0];
        }
    }
}

double randNorm(void) {
    // Implements the Box-Muller transform to sample from standard normal distribution
    double a, b, c;
    while (1) {
        a = 2.0*(double)rand()/RAND_MAX - 1.0;
        b = 2.0*(double)rand()/RAND_MAX - 1.0;
        c = a*a + b*b;
        if (c < 1.0) break;
    }
    return a*sqrt((-2.0*log(c))/c);
}


/* 
 * Function: scatter
 * ----------------
 * Multiple Coulomb scattering model from J.D. Jackson "Classical Electrodynamics" 3rd Ed. 
 * Modifies velocity components in place, preserving magnitude but changing direction.
 *
 *  vx: pointer to x velocity of proton
 *  vy: pointer to y velocity of proton
 *  vz: pointer to z velocity of proton
 *  Z: atomic number of the scattering material
 *  N: number density of the scattering material in cm^-3
 *  ds: step size through the material in cm
 */

void scatter(double *vx, double *vy, double *vz, double Z, double N, double ds) {
    double phi = 2*M_PI*(double)rand()/RAND_MAX; // Azimuthal angle (using Jackson's notation)
    double v = sqrt(pow(*vx,2)+pow(*vy,2)+pow(*vz,2));
    double p = (1.0/sqrt(1-pow(v/(C_LIGHT),2))) * (PROT_MASS) * v;
    double sqrtMeanSqAngle = 4*sqrt(M_PI*N*log(204*pow(Z,-1./3.))*ds)*Z*pow(PROT_CHARGE,2)/(p*v); // Eqn 13.66 from Jackson
    double theta = randNorm()*sqrtMeanSqAngle; // Polar (zenith) angle
    
    double sint = sin(theta);
    double cost = cos(theta);
    double sinp = sin(phi);
    double cosp = cos(phi);
    
    // Unit vector of the initial velocity
    double ux = *vx/v;
    double uy = *vy/v;
    double uz = *vz/v;

    // Compute orthonormal vector to initial velocity, w, by crossing u with a random vector r.
    double rx = (double)rand()/RAND_MAX;
    double ry = (double)rand()/RAND_MAX;
    double rz = (double)rand()/RAND_MAX;
    double wx = uy*rz-uz*ry;
    double wy = uz*rx-ux*rz;
    double wz = ux*ry-uy*rx;
    
    // In case r happened to be parallel to v, recompute w with new r. 
    // Super unlikely this ever actually loops.
    while (wx == 0.0 && wy == 0.0 && wz == 0.0) {
        rx = (double)rand()/RAND_MAX;
        ry = (double)rand()/RAND_MAX;
        rz = (double)rand()/RAND_MAX;
        wx = uy*rz-uz*ry;
        wy = uz*rx-ux*rz;
        wz = ux*ry-uy*rx;
    }

    double w_mag = sqrt(wx*wx+wy*wy+wz*wz);
    wx /= w_mag;
    wy /= w_mag;
    wz /= w_mag;

    // Multiply v by combined rotation matrix R = Rp*Rt, where Rt rotates 
    // angle theta about w, and Rp rotates angle phi about u.
    
    double vx_old = *vx;
    double vy_old = *vy;
    double vz_old = *vz;
    
    *vx = vx_old * ((cosp+(1-cosp)*ux*ux)*(cost+(1-cost)*wx*wx)
                 + ((1-cosp)*uy*ux-sinp*uz)*((1-cost)*wy*wx+sint*wz)
                 + ((1-cosp)*uz*ux+sinp*uy)*((1-cost)*wz*wx-sint*wy))
         +vy_old * ((cosp+(1-cosp)*ux*ux)*((1-cost)*wy*wx-sint*wz)
                 + ((1-cosp)*uy*ux-sinp*uz)*(cost+(1-cost)*wy*wy)
                 + ((1-cosp)*uz*ux+sinp*uy)*((1-cost)*wz*wy+sint*wx))
         +vz_old * ((cosp+(1-cosp)*ux*ux)*((1-cost)*wz*wx+sint*wy)
                 + ((1-cosp)*uy*ux-sinp*uz)*((1-cost)*wz*wy-sint*wx)
                 + ((1-cosp)*uz*ux+sinp*uy)*(cost+(1-cost)*wz*wz));

    *vy = vx_old * (((1-cosp)*uy*ux+sinp*uz)*(cost+(1-cost)*wx*wx)
                 + (cosp+(1-cosp)*uy*uy)*((1-cost)*wy*wx+sint*wz)
                 + ((1-cosp)*uz*uy-sinp*ux)*((1-cost)*wz*wx-sint*wy))
         +vy_old * (((1-cosp)*uy*ux+sinp*uz)*((1-cost)*wy*wx-sint*wz)
                 + (cosp+(1-cosp)*uy*uy)*(cost+(1-cost)*wy*wy)
                 + ((1-cosp)*uz*uy-sinp*ux)*((1-cost)*wz*wy+sint*wx))
         +vz_old * (((1-cosp)*uy*ux+sinp*uz)*((1-cost)*wz*wx+sint*wy)
                 + (cosp+(1-cosp)*uy*uy)*((1-cost)*wz*wy-sint*wx)
                 + ((1-cosp)*uz*uy-sinp*ux)*(cost+(1-cost)*wz*wz));

    *vz = vx_old * (((1-cosp)*uz*ux-sinp*uy)*(cost+(1-cost)*wx*wx)
                 + ((1-cosp)*uz*uy+sinp*ux)*((1-cost)*wy*wx+sint*wz)
                 + (cosp+(1-cosp)*uz*uz)*((1-cost)*wz*wx-sint*wy))
         +vy_old * (((1-cosp)*uz*ux-sinp*uy)*((1-cost)*wy*wx-sint*wz)
                 + ((1-cosp)*uz*uy+sinp*ux)*(cost+(1-cost)*wy*wy)
                 + (cosp+(1-cosp)*uz*uz)*((1-cost)*wz*wy+sint*wx))
         +vz_old * (((1-cosp)*uz*ux-sinp*uy)*((1-cost)*wz*wx+sint*wy)
                 + ((1-cosp)*uz*uy+sinp*ux)*((1-cost)*wz*wy-sint*wx)
                 + (cosp+(1-cosp)*uz*uz)*(cost+(1-cost)*wz*wz));
    
}

void pmover(int NP, int ns, double qm, double ds, double *x, double *y, double *z, 
		    double *vx, double *vy, double *vz, double *prop_dir, int ngridx, int ngridy, int ngridz,
	       	double *gridvals_flat, double *grid_offset, double *grid_spacing, 
		    int ntraces, double *traces_flat, int cyl_coords) {
	
    srand(time(NULL));


	double (*gridvals)[ngridy][ngridz][NUM_FIELDS] = (double (*)[ngridy][ngridz][NUM_FIELDS]) &gridvals_flat[0];
	double (*traces)[ns+3][3] = (double (*)[ns+3][3]) &traces_flat[0];
    
    #pragma omp parallel for
    for (int n = 0; n < NP; n++) {
		for (int step = 0; step < ns; step++) {

            if (n < ntraces) {
				traces[n][step+1][0] = x[n];
				traces[n][step+1][1] = y[n];
				traces[n][step+1][2] = z[n];
			}
            
            // Velocity along propogation axis
            double vs = vx[n]*prop_dir[0] + vy[n]*prop_dir[1] + vz[n]*prop_dir[2];
            
			x[n] += vx[n]*ds/vs;
			y[n] += vy[n]*ds/vs;
			z[n] += vz[n]*ds/vs;

			double* fieldvals = nn_interp(x[n], y[n], z[n], ngridx, ngridy, ngridz, gridvals, grid_offset, grid_spacing, cyl_coords);
            
            if (fieldvals[MATZ_IDX] > 0.0) {
                scatter(&vx[n], &vy[n], &vz[n], fieldvals[MATZ_IDX], fieldvals[DENS_IDX], ds);
            }
            
			vx[n] += (qm*ds/vs) * (fieldvals[EX_IDX] + (vy[n]*fieldvals[BZ_IDX] - vz[n]*fieldvals[BY_IDX])/C_LIGHT);
			vy[n] += (qm*ds/vs) * (fieldvals[EY_IDX] + (vz[n]*fieldvals[BX_IDX] - vx[n]*fieldvals[BZ_IDX])/C_LIGHT);
			vz[n] += (qm*ds/vs) * (fieldvals[EZ_IDX] + (vx[n]*fieldvals[BY_IDX] - vy[n]*fieldvals[BX_IDX])/C_LIGHT);

		}
	}
}


