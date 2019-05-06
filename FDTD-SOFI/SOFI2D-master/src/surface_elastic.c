/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 * 
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or 
  * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
 *   stress free surface condition
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void surface_elastic(int ndepth, int * gx, float ** vx, float ** vy, float ** sxx, float ** syy,
float **sxy, float  **pi, float  **u, float *hc, float * K_x, float * a_x, float * b_x, float ** psi_vxx){


	int i,j,m,fdoh,h1;
	float fjm, g, dhi;
	float  vxx, vyy;
	extern float DT, DH;
	extern int NX, FDORDER;
	extern int FW, BOUNDARY;
        extern int NPROCX, POS[3]; 
	extern int ABS_TYPE;
	fdoh = FDORDER/2;
	dhi = 1.0/DH;


	j=ndepth;     /* The free surface is located exactly in y=dh !! */
	for (i=1;i<=NX;i++){
		
		/*Mirroring the components of the stress tensor to make
			a stress free surface (method of imaging)*/
		syy[j][i]=0.0;

		
		vxx = 0.0;
		vyy = 0.0;
		 
		for (m=1; m<=fdoh; m++) {
		
			/*Mirroring the components of the stress tensor to make
			a stress free surface (method of imaging)*/
			syy[j-m][i]=-syy[j+m][i];
			sxy[j-m][i]=-sxy[j+m-1][i];
		       
			vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]);
			vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]);
		}
		 
		vxx *= dhi;
		vyy *= dhi;
	
		if (ABS_TYPE==1){
	      	/* apply PML boundary */    
             		/* left boundary */
             		if((!BOUNDARY) && (POS[1]==0) && (i<=FW)){
                        
                        	psi_vxx[j][i] = b_x[i] * psi_vxx[j][i] + a_x[i] * vxx;
                        	vxx = vxx / K_x[i] + psi_vxx[j][i];                 
             		}

             	/* right boundary */
             		if((!BOUNDARY) && (POS[1]==NPROCX-1) && (i>=NX-FW+1)){
                
                        	h1 = (i-NX+2*FW);
                        
                        	psi_vxx[j][h1] = b_x[h1] * psi_vxx[j][h1] + a_x[h1] * vxx;
                        	vxx = vxx / K_x[h1] + psi_vxx[j][h1];                                            
             		}  
             	}

		fjm=u[j][i]*2.0;
		g=pi[j][i];
		
		 /*Update sxx without vertical derivates (last update will be canceld)
                 *sxx=sxx_new - sxx =  DT*fjm*(2-fjm/g)vxx   -  ( g* ( vxx+vyy ) - fjm *vyy)
                 *                  = -(DT*(g-fmj)*(g-fmj)*vxx/g)-(DT*(g-fjm)*vyy) */

		sxx[j][i]+= -(DT*(g-fjm)*(g-fjm)*vxx/g)-(DT*(g-fjm)*vyy);
		
	}
	

	
}
