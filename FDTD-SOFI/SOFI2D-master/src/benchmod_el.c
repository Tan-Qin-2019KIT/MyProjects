/* -------------------------------------------------------------
 *   Model for overnight built, fullspace, elastic case
 *
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_elastic(float  **  rho, float **  pi, float **  u){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern int NX, NY, NXG, NYG,  POS[3], MYID;
	extern int WRITE_MODELFILES;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float muv, piv, vp, vs, rhov;
	int i, j, ii, jj;
	char modfile[STRING_SIZE];


	/*-----------------material property definition -------------------------*/	

	/* parameters for layer 1 */
	const float vp1=5100.0, vs1=2800.0, rho1=2100.0;


	/*-----------------------------------------------------------------------*/

	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (j=1;j<=NYG;j++){

			vp=vp1; vs=vs1; rhov=rho1; 

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

	/* each PE writes his model to disk */

	/* only the density model is written to file */
	if (WRITE_MODELFILES==2) {
		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}

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

		sprintf(modfile,"%s.SOFI2D.rho",MFILE);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MYID==0) mergemod(modfile,3);
	}
}



