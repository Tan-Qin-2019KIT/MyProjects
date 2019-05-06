/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
* Define damping profiles for CPML boundary condition
* This C-PML implementation is adapted from the 2nd order isotropic CPML code
* by Dimitri Komatitsch and based in part on formulas given in Roden and Gedney (2000). 
* 
* References:
* @ARTICLE{KoMa07,
* author = {Dimitri Komatitsch and Roland Martin},
* title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
*          at grazing incidence for the seismic wave equation},
* journal = {Geophysics},
* year = {2007},
* volume = {72},
* number = {5},
* pages = {SM155-SM167},
* doi = {10.1190/1.2757586}}
*
* @ARTICLE{MaKoEz08,
* author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
* title = {An unsplit convolutional perfectly matched layer improved at grazing
* incidence for seismic wave equation in poroelastic media},
* journal = {Geophysics},
* year = {2008},
* volume = {73},
* pages = {T51-T61},
* number = {4},
* doi = {10.1190/1.2939484}}
*
* @ARTICLE{MaKo09,
* author = {Roland Martin and Dimitri Komatitsch},
* title = {An unsplit convolutional perfectly matched layer technique improved
* at grazing incidence for the viscoelastic wave equation},
* journal = {Geophysical Journal International},
* year = {2009},
* volume = {179},
* pages = {333-344},
* number = {1},
* doi = {10.1111/j.1365-246X.2009.04278.x}}
*
* @ARTICLE{MaKoGe08,
* author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
* title = {A variational formulation of a stabilized unsplit convolutional perfectly
* matched layer for the isotropic or anisotropic seismic wave equation},
* journal = {Computer Modeling in Engineering and Sciences},
* year = {2008},
* volume = {37},
* pages = {274-304},
* number = {3}}
*
* The original CPML technique for Maxwell's equations is described in:
*
* @ARTICLE{RoGe00,
* author = {J. A. Roden and S. D. Gedney},
* title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
*          of the {CFS}-{PML} for Arbitrary Media},
* journal = {Microwave and Optical Technology Letters},
* year = {2000},
* volume = {27},
* number = {5},
* pages = {334-339},
* doi = {10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A}}
* extended version of 2D fdveps code
* S. Butzer 2010
*/

#include "fd.h"

void CPML_coeff(float * K_x, float * alpha_prime_x, float * a_x, float * b_x, 
            float * K_x_half, float * alpha_prime_x_half, float * a_x_half, float * b_x_half,
            float * K_y, float * alpha_prime_y, float * a_y, float * b_y, 
            float * K_y_half, float * alpha_prime_y_half, float * a_y_half, float * b_y_half,
            float * K_z, float * alpha_prime_z, float * a_z, float * b_z, 
            float * K_z_half, float * alpha_prime_z_half, float * a_z_half, float * b_z_half)
{

	extern float DX, DY, DZ , VPPML, DT, FPML;
	extern int FW;
	extern FILE *FP;
	
	int i, l,i1;
	      
      float * d_x=NULL,* d_x_half=NULL,* d_y=NULL,* d_y_half=NULL,* d_z=NULL,* d_z_half=NULL;
      
      const float npower = 3;  /*  power to compute d0 profile */
      const float k_max_PML = 2.0;   /* (from Gedney page 8.11) */
      const float alpha_max_PML = 2.0 * PI * (FPML/2.0);   /* from festa and Vilotte 2.0*...*/
      const float Rcoef = 0.0008;       /* reflection coefficient (INRIA report section 6.1) */
      const float a = 0.25, b = 0.75 , c = 0.0; 
      float d0_x, d0_y, d0_z, position_norm, position_in_PML;

      d_x = vector(1,2*FW);
      d_x_half = vector(1,2*FW);
      d_y = vector(1,2*FW);
      d_y_half = vector(1,2*FW);
      d_z = vector(1,2*FW);
      d_z_half = vector(1,2*FW);
     
      /* compute d0 from INRIA report section 6.1 */
      d0_x = - (npower + 1) * VPPML * log(Rcoef) / (2.0 * FW*DX);
      d0_y = - (npower + 1) * VPPML * log(Rcoef) / (2.0 * FW*DY);
      d0_z = - (npower + 1) * VPPML * log(Rcoef) / (2.0 * FW*DZ);


	
      /* damping in the X direction */
      /* -------------------------- */

   	for (i=1;i<=FW+1;i++){
	
        K_x[i] = 1.0;

        /* define damping profile at the grid points */
        position_in_PML = (FW-i+1)*DX; /*distance to boundary in meter */
        position_norm = position_in_PML / (FW*DX); /*normalised by PML thickness*/

        d_x[i] = d0_x *(a*position_norm+b*pow(position_norm,npower)+c*pow(position_norm,(4)));

        /* this taken from Gedney page 8.2 */
        K_x[i] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
        alpha_prime_x[i] = alpha_max_PML * (1.0 - position_norm);

	if(alpha_prime_x[i] < 0.0){ fprintf(FP,"ERROR:alpha_prime_x[i] < 0.0, i %d", i);}
	
	b_x[i] = exp(- (d_x[i] / K_x[i] + alpha_prime_x[i]) * DT);

 	/* avoid division by zero outside the PML */
        if(abs(d_x[i]) > 1.0e-6){ a_x[i] = d_x[i] * (b_x[i] - 1.0) / (K_x[i] * (d_x[i] + K_x[i] * alpha_prime_x[i]));}
	else a_x[i]=0.0;
	
	if(i<=FW){

        /* define damping profile at half the grid points (half a grid point in -x)*/
        position_in_PML = (FW-i+0.5) *DX;
        position_norm = position_in_PML / (FW*DX);

	i1=i;
	
	K_x_half[i1] = 1.0;
        d_x_half[i1] = d0_x * (a*position_norm+b*pow(position_norm,npower)+c*pow(position_norm,(4)));

        if(position_in_PML < 0.0) {fprintf(FP,"ERROR: Position in PML (x-boundary) smaller 0");}
          
        /* this taken from Gedney page 8.2 */
        K_x_half[i1] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
        alpha_prime_x_half[i1] = alpha_max_PML * (1.0 - position_norm);
           
        /* just in case, for -5 at the end */
        if(alpha_prime_x_half[i1] < 0.0) {fprintf(FP,"ERROR:alpha_prime_x_half[i] < 0.0, i %d", i);}

        b_x_half[i1] = exp(- (d_x_half[i1] / K_x_half[i1] + alpha_prime_x_half[i1]) * DT);

        if(abs(d_x_half[i1]) > 1.0e-6){ a_x_half[i1] = d_x_half[i1] * (b_x_half[i1] - 1.0) / (K_x_half[i1] * (d_x_half[i1] + K_x_half[i1] * alpha_prime_x_half[i1]));}

	/* right boundary --> mirroring left boundary*/
	
	l = 2* FW -i+1;

	if(i>1){
	K_x[l+1]=K_x[i];
	b_x[l+1] = b_x[i];
	if(abs(d_x[i]) > 1.0e-6){ a_x[l+1] = a_x[i]; }
	}

	K_x_half[l]=K_x_half[i];
        b_x_half[l] = b_x_half[i];  /*half a grid point in +x)*/
        if(abs(d_x[i]) > 1.0e-6){ a_x_half[l] = a_x_half[i]; }

        } 
	}



      /* damping in the Y direction */
      /* -------------------------- */

        for (i=1;i<=FW+1;i++){
	
        K_y[i] = 1.0; 
          
        /* define damping profile at the grid points */
        position_in_PML = (FW-i+1)*DY; /*distance to boundary in meter */
        position_norm = position_in_PML / (FW*DY); /*normalised by PML thickness*/

        d_y[i] = d0_y * (a*position_norm+b*pow(position_norm,npower)+c*pow(position_norm,(4)));

        /* this taken from Gedney page 8.2 */
        K_y[i] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
        alpha_prime_y[i] = alpha_max_PML * (1.0 - position_norm);

	/* just in case, for -5 at the end */
        if(alpha_prime_y[i] < 0.0){ fprintf(FP,"ERROR:alpha_prime_y[i] < 0.0, i %d", i);}

	b_y[i] = exp(- (d_y[i] / K_y[i] + alpha_prime_y[i]) * DT);

 	/* avoid division by zero outside the PML */
        if(abs(d_y[i]) > 1.0e-6){ a_y[i] = d_y[i] * (b_y[i] - 1.0) / (K_y[i] * (d_y[i] + K_y[i] * alpha_prime_y[i]));}
      	else a_x[i]=0.0;

	if(i<=FW){

          /* define damping profile at half the grid points (half a grid point in -x)*/
        position_in_PML = (FW-i+0.5) *DY;
        position_norm = position_in_PML / (FW*DY);

	i1=i;
	K_y_half[i1] = 1.0;
        d_y_half[i1] = d0_y * (a*position_norm+b*pow(position_norm,npower)+c*pow(position_norm,(4)));

        if(position_in_PML < 0.0) {fprintf(FP,"ERROR: Position in PML (y-boundary) <0");}
          
        /* this taken from Gedney page 8.2 */
        K_y_half[i1] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
        alpha_prime_y_half[i1] = alpha_max_PML * (1.0 - position_norm);
      
        if(alpha_prime_y_half[i1] < 0.0) {fprintf(FP,"ERROR:alpha_prime_y_half[i] < 0.0, i %d", i);}
        b_y_half[i1] = exp(- (d_y_half[i1] / K_y_half[i1] + alpha_prime_y_half[i1]) * DT);
          
      	if(abs(d_y_half[i1]) > 1.0e-6){ a_y_half[i1] = d_y_half[i1] * (b_y_half[i1] - 1.0) / (K_y_half[i1] * (d_y_half[i1] + K_y_half[i1] * alpha_prime_y_half[i1]));}
	
        /* right boundary --> mirroring left boundary*/
        l = 2* FW -i+1;
	
	if(i>1){
	K_y[l+1] = K_y[i];
	b_y[l+1] = b_y[i];
	if(abs(d_y[i]) > 1.0e-6){ a_y[l+1] = a_y[i]; }
	}
  
	
	K_y_half[l]=K_y_half[i];
        b_y_half[l] = b_y_half[i];  /*half a grid point in +x)*/ 
        if(abs(d_y[i]) > 1.0e-6){ a_y_half[l] = a_y_half[i]; }
	}
        } 


       /* damping in the Z direction */
      /* -------------------------- */

	for (i=1;i<=FW+1;i++){
	
        K_z[i] = 1.0;
                    
    	/* define damping profile at the grid points */
      	position_in_PML = (FW-i+1)*DZ; /*distance to boundary in meter */
      	position_norm = position_in_PML / (FW*DZ); /*normalised by PML thickness*/

      	d_z[i] = d0_z * (a*position_norm+b*pow(position_norm,npower)+c*pow(position_norm,(4)));

  	/* this taken from Gedney page 8.2 */
    	K_z[i] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
    	alpha_prime_z[i] = alpha_max_PML * (1.0 - position_norm);
           
	/* just in case, for -5 at the end */
    	if(alpha_prime_z[i] < 0.0){ fprintf(FP,"ERROR:alpha_prime_z[i] < 0.0, i %d", i);}

	b_z[i] = exp(- (d_z[i] / K_z[i] + alpha_prime_z[i]) * DT);

 	/* avoid division by zero outside the PML */
        if(abs(d_z[i]) > 1.0e-6){ a_z[i] = d_z[i] * (b_z[i] - 1.0) / (K_z[i] * (d_z[i] + K_z[i] * alpha_prime_z[i]));}
	
	if(i<=FW){
    	/* define damping profile at half the grid points (half a grid point in -x)*/
        position_in_PML = (FW-i+0.5) *DZ;
        position_norm = position_in_PML / (FW*DZ);

	i1=i;

	K_z_half[i1] = 1.0;
        d_z_half[i1] = d0_z *(a*position_norm+b* pow(position_norm,npower)+c*pow(position_norm,(4)));

        if(position_in_PML < 0.0) {fprintf(FP,"ERROR: Position in PML (y-boundary) <0");}
 
        /* this taken from Gedney page 8.2 */
        K_z_half[i1] = 1.0 + (k_max_PML - 1.0) * pow(position_norm,npower);
        alpha_prime_z_half[i1] = alpha_max_PML * (1.0 - position_norm);

        if(alpha_prime_z_half[i1] < 0.0) {fprintf(FP,"ERROR:alpha_prime_z_half[i] < 0.0, i %d", i);}
 
        b_z_half[i1] = exp(- (d_z_half[i1] / K_z_half[i1] + alpha_prime_z_half[i1]) * DT);

        if(abs(d_z_half[i1]) > 1.0e-6){ a_z_half[i1] = d_z_half[i1] * (b_z_half[i1] - 1.0) / (K_z_half[i1] * (d_z_half[i1] + K_z_half[i1] * alpha_prime_z_half[i1]));}
	
        /* right boundary --> mirroring left boundary*/
        l = 2* FW -i+1;

	if(i>1){
	K_z[l+1] = K_z[i];
	b_z[l+1] = b_z[i];
	if(abs(d_z[i]) > 1.0e-6){ a_z[l+1] = a_z[i]; }
	}
	
	
	K_z_half[l]=K_z_half[i];
        b_z_half[l] = b_z_half[i];  /*half a grid point in +x)*/
        if(abs(d_z[i]) > 1.0e-6){ a_z_half[l] = a_z_half[i]; }
	}
}
    free_vector(d_x,1,2*FW);
    free_vector(d_x_half,1,2*FW);
    free_vector(d_y,1,2*FW);
    free_vector(d_y_half,1,2*FW);
    free_vector(d_z,1,2*FW);
    free_vector(d_z_half,1,2*FW);
}



