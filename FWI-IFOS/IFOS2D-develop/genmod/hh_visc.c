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

/*
 *   Model layer over half space (homogeneous)
 *   update 11.04.02, T. Bohlen
 *   last update 11. Okt. 2011 M. Schaefer
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, 
float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU, DH;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	

		/* local variables */
	float rhov, vp, vs, y;
	float *pts;
	int i, j, l, ii, jj;
	 
	
	/* parameters for layer 1 */
	const float vp1=450.0, vs1=270.0, rho1=1800.0, h=5.0;
	
	/* parameters for layer 2 */
	const float vp2=450.0, vs2=270.0, rho2=1800.0;
	
	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

			

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				y=(float)j*DH;
			
				if (y<=h){
				 vp=vp1; vs=vs1; rhov=rho1; }

				
				 else{
 				 vp=vp2; vs=vs2; rhov=rho2; }
                    
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					taus[jj][ii]=TAU;
					taup[jj][ii]=TAU;
					u[jj][ii]=vs;
					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
				}
			}
		}	

		

	
	/* each PE writes his model to disk */

	writemod(MFILE,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(MFILE,3);

	free_vector(pts,1,L);
}



