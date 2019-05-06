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

/* $Id: hh_el.c,v 2.3 2007/08/21 13:16:19 tbohlen Exp $*/
/*
 *   Model defined by flnode file. In this version a constant Q-model is used. It is defined via the 
 *   relaxation frequencies and the TAU value given in the input file. Nevertheless the flnode file
 *   must contain two comlumns with tau values for taup and taus. 
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], MYID, L;
	extern float DH, DT, *FL, TAU;
	extern char  MFILE[STRING_SIZE];	

		/* local variables */
	float muv, piv, vp, vs, rhov, ts, tp, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE];
	
	FILE *flfile;
	int nodes;
	char cline[256];
	
	float *fldepth, *flrho, *flvp, *flvs, *fltaus, *fltaup;
	
	nodes=7;
	
	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);
	flvs=vector(1,nodes);
	fltaup=vector(1,nodes);
	fltaus=vector(1,nodes);
	
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}
	 
	
	
	/*read FL nodes from File*/
	
	flfile=fopen("model/final.mod.flnodes.Q20","r");
	if (flfile==NULL) declare_error(" FL-file could not be opened !");
	
	
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &flvs[l], &fltaup[l], &fltaus[l]);
		}
		else l=l-1;
	
	}
	
	
		
	/*-----------------------------------------------------------------------*/

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (l=1;l<nodes;l++){
				if (fldepth[l]==fldepth[l+1]){
					if ((i==1) && (MYID==0)){
					printf("depth: %f m: double node\n",fldepth[l]);}}
				else{
					for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
			
						vp=0.0; vs=0.0; rhov=0.0;
				  
						vp=(DH*(j-1)-fldepth[l])*(flvp[l+1]-flvp[l])/(fldepth[l+1]-fldepth[l])+flvp[l];
						vp=vp*1000.0;
						vs=(DH*(j-1)-fldepth[l])*(flvs[l+1]-flvs[l])/(fldepth[l+1]-fldepth[l])+flvs[l];
						vs=vs*1000.0;
						rhov=(DH*(j-1)-fldepth[l])*(flrho[l+1]-flrho[l])/(fldepth[l+1]-fldepth[l])+flrho[l];
						rhov=rhov*1000.0;
						
						muv=vs;
						piv=vp;
						ts=TAU;
						tp=TAU;
						
						/* only the PE which belongs to the current global gridpoint 
				  		is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
				    		(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						u[jj][ii]=muv;
						rho[jj][ii]=rhov;
						pi[jj][ii]=piv;
						taus[jj][ii]=ts;
						taup[jj][ii]=tp;
						}
					}
				}
			}
			
			for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			  
				vp=0.0; vs=0.0; rhov=0.0;
				vp=flvp[nodes]*1000.0; vs=flvs[nodes]*1000.0; rhov=flrho[nodes]*1000.0;
				
				muv=vs;
				piv=vp;
				ts=TAU;
				tp=TAU;

				/* only the PE which belongs to the current global gridpoint 
				  		is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
				    		(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						u[jj][ii]=muv;
						rho[jj][ii]=rhov;
						pi[jj][ii]=piv;
						taus[jj][ii]=ts;
						taup[jj][ii]=tp;
						}
			}
		}
		
		
        /* each PE writes his model to disk */

	/* each PE writes his model to disk */
	sprintf(modfile,"%s.mu",MFILE);
	writemod(modfile,u,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);	
	
	
	/* each PE writes his model to disk */
	sprintf(modfile,"%s.pi",MFILE);
	writemod(modfile,pi,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);	
	
	
	
	/* each PE writes his model to disk */
	sprintf(modfile,"%s.rho",MFILE);
	writemod(modfile,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
	
	
	
	/* each PE writes his model to disk */
	sprintf(modfile,"%s.taus",MFILE);
	writemod(modfile,taus,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
	
	
	
	/* each PE writes his model to disk */
	sprintf(modfile,"%s.taup",MFILE);
	writemod(modfile,taup,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(modfile,3);
	
	free_vector(fldepth,1,nodes);
	free_vector(flrho,1,nodes);
	free_vector(flvp,1,nodes);
	free_vector(flvs,1,nodes);
	free_vector(fltaup,1,nodes);
	free_vector(fltaus,1,nodes);
	
	free_vector(pts,1,L);
	 
}



