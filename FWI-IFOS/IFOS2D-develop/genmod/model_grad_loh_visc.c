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

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern char INV_MODELFILE[STRING_SIZE];
	extern float DH, *FL, TAU, DT;
		/* local variables */
	float vp, vs, rhov, grad, y1, y2, ts, tp, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE]; 
	
	/* gradient for Vs: parameters for velocity in m/s, depth h in meter */
	const float vs1=120.0, vs2=300.0, h=12.0;  
	 
	/* layer over halfspace for Vp and rho: parameters parameters for velocity in m/s, depth of layer in meter*/
	const float vp1=300.0, vp2=1800.0, rho1=1700.0, rho2=2000.0, layer=6.0;
	
	/*-----------------------------------------------------------------------*/
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
	        eta[l]=DT/pts[l];
	}
	
	
	y1=h/DH;
	y2=layer/DH;
	if(y1==NYG) declare_error(" \n y is equal NYG !! see src/model_grad.c  \n ");
	grad=(vs2-vs1)/y1;	
	
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			
				/* gradient for Vs */
				if(j<=y1){
				vs=vs1+(j*grad);
				}
				
				else{				
				vs=vs2;
				}
				
				/* layer over halfspace for Vp and rho*/
				if(j<=y2){
				vp=vp1;
				rhov=rho1;
				}
				
				else{
				vp=vp2;
				rhov=rho2;
				}
				
				ts=TAU;
				tp=TAU;
				
				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=vs;
					rho[jj][ii]=rhov;
					pi[jj][ii]=vp;
					taus[jj][ii]=ts;
					taup[jj][ii]=tp;
				}
			}
		}	

		
sprintf(modfile,"%s_rho_it_0.bin",INV_MODELFILE);
writemod(modfile,rho,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_0.bin",INV_MODELFILE);
writemod(modfile,u,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vp_it_0.bin",INV_MODELFILE);
writemod(modfile,pi,3);
MPI_Barrier(MPI_COMM_WORLD);
if (MYID==0) mergemod(modfile,3);

free_vector(pts,1,L);
}



