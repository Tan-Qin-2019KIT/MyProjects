/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2005  <name of author>
 *
 * This file is part of DENISE.
 * 
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Read elastic model properties (vp,vs,density) from files  
 *
 *  Copyright (c)  T. Bohlen
 *  last update 29.06.2003
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
   to read data from model-files for viscoelastic simulation */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	extern float DH;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, piv, vp, vs, x, y;
	float vp0, vs0, rho0, vp1, vs1, rho1, r, r1;
	int i, j, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_rho;
	char filename[STRING_SIZE];
        
	 /* location of the spheres */
        const float X01 = 80.0;
        const float Y01 = 130.0;

	
	vp0=2000.0;
	vs0=0.0;
	rho0=1000.0;
	
	vp1=1700.0;
	vs1=0.0;
	rho1=1000.0;
	
	r1=10.0;
	

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
			
			   x=(float)i*DH;
                           y=(float)j*DH;

			   vp = vp0;
			   vs = vs0;
			   rhov = rho0;
			
			   r=sqrt(((x-X01)*(x-X01))+((y-Y01)*(y-Y01)));
                                
			   if(r<=r1){
                              vp=vp1; vs=vs1; rhov=rho1;
                           }

				
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
                                
				u[jj][ii]=vs;
                                rho[jj][ii]=rhov;
                                pi[jj][ii]=vp;
				
				}
			}
		}
	

	
	/* each PE writes his model to disk */
	   
	   
	sprintf(filename,"%s.fdveps.vp",MFILE);

	writemod(filename,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
	
	
	sprintf(filename,"%s.fdveps.vs",MFILE);

	writemod(filename,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);

        sprintf(filename,"%s.fdveps.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);
}




