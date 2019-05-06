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
/* -------------------------------------------------------------
 *   Model homogeneous half space
 *
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_visco(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU, DH;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern int WRITE_MODELFILES;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float Rhov, muv, piv, Vp, Vs, y;
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, ii, jj;
	char modfile[STRING_SIZE];	
	float ** pwavemod=NULL, ** swavemod=NULL;

	/*-----------------material property definition -------------------------*/	

	/* parameters for layer 1 */
	const float vp1=3500.0, vs1=2000.0, rho1=2000.0, h=10000.0;

	/* parameters for layer 2 */
	const float vp2=5400.0, vs2=3700.0, rho2=2500.0;


	/*-----------------------------------------------------------------------*/

	if (WRITE_MODELFILES==1) {
		pwavemod  =  matrix(0,NY+1,0,NX+1);
		swavemod  =  matrix(0,NY+1,0,NX+1);
	}

	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	ws=2.0*PI*FL[1];

	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}



	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){

			/* calculate coordinate in m */
			y=(float)j*DH;

			/* two layer case */
			if (y<=h){
				Vp=vp1; Vs=vs1; Rhov=rho1; }


			else{
				Vp=vp2; Vs=vs2; Rhov=rho2; }


			muv=Vs*Vs*Rhov/(1.0+sumu);
			piv=Vp*Vp*Rhov/(1.0+sumpi);

			/* only the PE which belongs to the current global gridpoint
				  is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) &&
					(POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				taus[jj][ii]=ts;
				taup[jj][ii]=tp;
				u[jj][ii]=muv;
				rho[jj][ii]=Rhov;
				pi[jj][ii]=piv;
				if (WRITE_MODELFILES==1)
				{
					pwavemod[jj][ii]=Vp;
					swavemod[jj][ii]=Vs;
				}
			}
		}
	}




	/* each PE writes his model to disk */

	/* all models are written to file */
	if (WRITE_MODELFILES==1) {
		sprintf(modfile,"%s.SOFI2D.u",MFILE);
		writemod(modfile,u,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.pi",MFILE);
		writemod(modfile,pi,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.ts",MFILE);
		writemod(modfile,taus,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.tp",MFILE);
		writemod(modfile,taup,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.vp",MFILE);
		writemod(modfile,pwavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.vs",MFILE);
		writemod(modfile,swavemod,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);

		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	/* only density is written to file */
	if (WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

	free_vector(pts,1,L);
	if (WRITE_MODELFILES==1) {
		free_matrix(pwavemod,0,NY+1,0,NX+1);
		free_matrix(swavemod,0,NY+1,0,NX+1);
	}
}
