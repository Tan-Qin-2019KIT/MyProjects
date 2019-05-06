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
 *   Read elastic model properties (vp,vs,density) from files
 *   This file contains function readmod, which has the purpose
 *   to read data from model-files for viscoelastic simulation
 *
 *  ----------------------------------------------------------------------*/


#include "fd.h"

void readmod_elastic(float  **  rho, float **  pi, float **  u){

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern char  MFILE[STRING_SIZE];	
	extern FILE *FP;

		
	/* local variables */
	float rhov, muv, piv, vp, vs;
	int i, j, ii, jj;
	FILE *fp_vs, *fp_vp, *fp_rho;
	char filename[STRING_SIZE];




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

	   

	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (j=1;j<=NYG;j++){
			fread(&vp, sizeof(float), 1, fp_vp);
			fread(&vs, sizeof(float), 1, fp_vs);
			fread(&rhov, sizeof(float), 1, fp_rho);
				
		if ((isnan(vp)) && (MYID==0)) {
           	declare_error(" Found NaN-Values in Vp-Model !");}

		if ((isnan(vs)) && (MYID==0)) {
                declare_error(" Found NaN-Values in Vs-Model !");}

		if ((isnan(rhov)) && (MYID==0)) {
                declare_error(" Found NaN-Values in Rho-Model !");}


			muv=vs*vs*rhov;
			piv=vp*vp*rhov;

			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;

					u[jj][ii]=muv;
					rho[jj][ii]=rhov;
					pi[jj][ii]=piv;
				}
			}
		}
	




	fclose(fp_vp);
	fclose(fp_vs);
	fclose(fp_rho);
	
	
	/* each PE writes his model to disk */
	   
	   
	sprintf(filename,"%s.sofi2D.rho",MFILE);

	writemod(filename,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(filename,3);

}




