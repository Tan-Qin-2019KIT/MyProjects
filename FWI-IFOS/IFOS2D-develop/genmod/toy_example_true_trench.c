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
 *   Model defined by flnode file. In this version a constant Q-model with Q_p=Q_s is used. It 
 *   is defined via the relaxation frequencies and the TAU value given in the input file.  
 */

#include "fd.h"

void model(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], MYID, L;
	extern float DH, DT, *FL, TAU;
	extern char  MFILE[STRING_SIZE];
	extern char INV_MODELFILE[STRING_SIZE];	

		/* local variables */
	float muv, piv, vp, vs, rhov, ts, tp, *pts;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE];
	
	FILE *flfile;
	int nodes;
	char cline[256];
	
	float *fldepth, *flrho, *flvp, *flvs;
	
	/*---------------------------------------
	 Trench "Ettlinger Linie" units in m  
	----------------------------------------*/

	/* trench=1 - elliptic shaped trench , trench=2 - V shaped trench */
	int trench=2;

	float width=10;
	float depth=3;	
	float pos=NXG/2*DH;

	// Parameters of trench "Ettlinger Linie" velocities in km/s and density in kg/m^3
	float evs=0.0917;
	float evp=0.2569;
	float erho=1.5;	
	
	/* fopen FL-node file*/
	flfile=fopen("model_true/flnodes.toy_example","r");
        if (flfile==NULL) declare_error(" FL-file could not be opened !");
	/* nodes specified in FL-node File*/
	nodes=7;

	/*-----------No changes below that line needed ---------------------------*/

	/* 2 linear lines as border for the v-shape */
	float m1=2*depth/width;
	float m2=-m1;
	float c1=depth*(1-(2*pos/width));
	float c2=depth*(1+(2*pos/width));

	

	
	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);
	flvs=vector(1,nodes);
		
	
	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}
	 
	
	
	/*read FL nodes from File*/
	
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &flvs[l]);
		}
		else l=l-1;
	
	}
	
	
	if(MYID==0){
	printf(" ------------------------------------------------------------------ \n\n");
	printf(" Information of FL nodes: \n\n");
	printf(" \t depth \t rho \t vp \t vs \n\n");
	
	for (l=1;l<=nodes;l++){
	printf(" \t %f \t %f \t %f \t %f\n\n",fldepth[l],flrho[l],flvp[l],flvs[l]);
	}
	printf(" ------------------------------------------------------------------ \n\n");
	}
	/*-----------------------------------------------------------------------*/
	
	
		
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
			/*Trench Ettlinger Linie*/
			if (trench==1) {
			 	for (j=1;j<=(int) depth/DH+1;j++){
				
					if ( ( ((pos-i*DH) * (pos-i*DH))/(width/2*width/2) + ((j*DH)*(j*DH))/(depth*depth) ) <=1){
                                	vp=evp*1000.0; 
					vs=evs*1000.0; 
					rhov=erho*1000.0;
				
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
			 /*Trench V-shape Ettlinger Linie*/
                        if (trench==2) {
                                for (j=1;j<=(int) depth/DH+1;j++){
					
                                        if ( ((i*DH)*m1+c1 >=0) && j*DH <= (i*DH)*m1+c1  && j*DH <= (i*DH)*m2+c2){
                                        vp=evp*1000.0;
                                        vs=evs*1000.0;
                                        rhov=erho*1000.0;

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

		}
		
		
        sprintf(modfile,"%s_rho_it0.bin",INV_MODELFILE);
	writemod(modfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_rho_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s_vs_it0.bin",INV_MODELFILE);
	writemod(modfile,u,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vs_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s_vp_it0.bin",INV_MODELFILE);
	writemod(modfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vp_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s_taus_it0.bin",INV_MODELFILE);
	writemod(modfile,taus,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_taus_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s_taup_it0.bin",INV_MODELFILE);
	writemod(modfile,taup,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_taup_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	
	free_vector(fldepth,1,nodes);
	free_vector(flrho,1,nodes);
	free_vector(flvp,1,nodes);
	free_vector(flvs,1,nodes);
	
	free_vector(pts,1,L);
	 
}



