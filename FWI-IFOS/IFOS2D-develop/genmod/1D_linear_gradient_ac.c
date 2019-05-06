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
 *   Model homogeneous half space
 *   last update 11.04.02, T. Bohlen
 */

#include "fd.h"

void model_acoustic(float  **  rho, float **  pi){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH;
		/* local variables */
	float vp, rhov, grad1, grad3, y;
	int i, j, ii, jj;
	char modfile[STRING_SIZE]; 
	
	/* parameters for layer 1 */
	const float vp1=500.0, rho1=1800.0, h=15.0;
	
	/* parameters for layer 2 due to calculation of grad1, grad2 and grad3*/
	const float vp2=1200.0, rho2=2000.0;
	
	/*-----------------------------------------------------------------------*/

	y=h/DH;
	if(y==NYG) declare_error(" \n y is equal NYG !! see src/model_grad.c  \n ");
	grad1=(vp2-vp1)/y;
	grad3=(rho2-rho1)/y;	
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				if(j<=y){
				vp=vp1+(j*grad1);
				rhov=rho1+(j*grad3);
				}
				
				else{				
				vp=vp2;
				rhov=rho2;
				}
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
				}
			}
		}	

		
sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
writemod(modfile,rho,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vp_it_0.bin",INV_MODELFILE);
writemod(modfile,pi,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);
}



