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
 * Read elastic model properties (vp,vs,density) from files
 * This file contains function readmod, which has the purpose
 * to read data from model-files for viscoelastic simulation
 *
 * ----------------------------------------------------------------------*/

#include "fd.h"

#include "fd.h"

void readmod_visco(float  **  rho, float **  pi, float **  u,
float **  taus, float **  taup, float *  eta){

	extern float DT, *FL, TAU, TS;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, piv, vp, vs, qp, qs;
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_rho, *fp_qp ,*fp_qs;
	char filename[STRING_SIZE];



	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}

	ts=TAU;  
	tp=TAU;

	//previously : ws=2.0*PI*FL[1];
	ws=2.0*PI/TS;
	
	sumu=0.0; 
	sumpi=0.0;
	for (l=1;l<=L;l++){
		sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
		sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
	}

		

	   fprintf(FP,"\n...reading model information from model-files...\n");

	   fprintf(FP,"\t P-wave velocities:\n\t %s.vp\n\n",MFILE);
	   sprintf(filename,"%s.vp",MFILE);
	   fp_vp=fopen(filename,"r");
	   if ((fp_vp==NULL) && (MYID==0)) declare_error(" Could not open model file for P velocities ! ");


	   fprintf(FP,"\t Shear wave velocities:\n\t %s.vs\n\n",MFILE);
	   sprintf(filename,"%s.vs",MFILE);
	   fp_vs=fopen(filename,"r");
	   if ((fp_vs==NULL) && (MYID==0)) declare_error(" Could not open model file for shear velocities ! ");

	   fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
	   sprintf(filename,"%s.rho",MFILE);
	   fp_rho=fopen(filename,"r");
	   if ((fp_rho==NULL) && (MYID==0)) declare_error(" Could not open model file for densities ! ");

	   fprintf(FP,"\t Qp:\n\t %s.qp\n\n",MFILE);
	   sprintf(filename,"%s.qp",MFILE);
	   fp_qp=fopen(filename,"r");
	   if ((fp_qp==NULL) && (MYID==0)) declare_error(" Could not open model file for Qp-values ! ");

	   fprintf(FP,"\t Qs:\n\t %s.qs\n\n",MFILE);
	   sprintf(filename,"%s.qs",MFILE);
	   fp_qs=fopen(filename,"r");
	   if ((fp_qs==NULL) && (MYID==0)) declare_error(" Could not open model file for Qs-values ! ");
	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
			fread(&qp, sizeof(float), 1, fp_qp);
			fread(&qs, sizeof(float), 1, fp_qs);
				
		if ((isnan(vp)) && (MYID==0)) {
                	declare_error(" Found NaN-Values in Vp-Model !");}

                if ((isnan(vs)) && (MYID==0)) {
                	declare_error(" Found NaN-Values in Vs-Model !");}

                if ((isnan(rhov)) && (MYID==0)) {
                	declare_error(" Found NaN-Values in Rho-Model !");}

		if ((isnan(qp)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in Vs-Model !");}

                if ((isnan(qs)) && (MYID==0)) {
                        declare_error(" Found NaN-Values in Rho-Model !");}


			muv=vs*vs*rhov/(1.0+sumu);
			piv=vp*vp*rhov/(1.0+sumpi);

			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					taus[jj][ii]=2.0/qs;
					taup[jj][ii]=2.0/qp;
					u[jj][ii]=muv;
					rho[jj][ii]=rhov;
					pi[jj][ii]=piv;
				}
			}
		}
	




	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
	fclose(fp_qp);
	fclose(fp_qs);
	
	
	/* each PE writes his model to disk */
	   
	   
	sprintf(filename,"%s.sofi2D.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);

	free_vector(pts,1,L);
}




