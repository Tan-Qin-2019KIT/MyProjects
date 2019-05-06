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

void model_acoustic(float **rho, float **pi){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG, POS[3], MYID;
	extern float DH, DT;
	extern char  MFILE[STRING_SIZE];
	extern char INV_MODELFILE[STRING_SIZE];	

		/* local variables */
	float piv, vp, rhov;
	int i, j, ii, jj, l;
	char modfile[STRING_SIZE];
	
	FILE *flfile;
	int nodes;
	char cline[256];
	
	float *fldepth, *flrho, *flvp;
	
	nodes=6;
	
	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);
	
	
	/* vector for maxwellbodies */
// 	pts=vector(1,L);
// 	for (l=1;l<=L;l++) {
// 		pts[l]=1.0/(2.0*PI*FL[l]);
// 		eta[l]=DT/pts[l];
// 	}	 
	
	
	/*read FL nodes from File*/
	
	flfile=fopen("model_true/flnodes.toy_example_ac","r");
	if (flfile==NULL) declare_error(" FL-file could not be opened !");
	
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f",&fldepth[l], &flrho[l], &flvp[l]);
		}
		else l=l-1;
	
	}
	
	
	if(MYID==0){
	printf(" ------------------------------------------------------------------ \n\n");
	printf(" Information of FL nodes: \n\n");
	printf(" \t depth \t rho \t vp \n\n");
	
	for (l=1;l<=nodes;l++){
	printf(" \t %f \t %f \t %f \n\n",fldepth[l],flrho[l],flvp[l]);
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
			
						vp=0.0; rhov=0.0;
				  
						vp=(DH*(j-1)-fldepth[l])*(flvp[l+1]-flvp[l])/(fldepth[l+1]-fldepth[l])+flvp[l];
						vp=vp*1000.0;
						rhov=(DH*(j-1)-fldepth[l])*(flrho[l+1]-flrho[l])/(fldepth[l+1]-fldepth[l])+flrho[l];
						rhov=rhov*1000.0;
						
						
						piv=vp;
						
						/* only the PE which belongs to the current global gridpoint 
				  		is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
				    		(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						rho[jj][ii]=rhov;
						pi[jj][ii]=piv;
						}
					}
				}
			}
			
			for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			  
				vp=0.0; rhov=0.0;
				vp=flvp[nodes]*1000.0; rhov=flrho[nodes]*1000.0;
				
				piv=vp;

				/* only the PE which belongs to the current global gridpoint 
				  		is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
				    		(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						rho[jj][ii]=rhov;
						pi[jj][ii]=piv;
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
	
	sprintf(modfile,"%s_vp_it0.bin",INV_MODELFILE);
	writemod(modfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vp_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	
	free_vector(fldepth,1,nodes);
	free_vector(flrho,1,nodes);
	free_vector(flvp,1,nodes);
	 
}
