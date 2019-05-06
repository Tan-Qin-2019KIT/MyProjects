/* -------------------------------------------------------------
 *   Model for overnight built, fullspace, viscoelastic case
 *
 *   ------------------------------------------------------------- */

#include "fd.h"

void model_visco(float  **  rho, float **  pi, float **  u, float **  taus, float **  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */

	extern float DT, *FL, TAU;
	extern int NX, NY, NXG, NYG,  POS[3], L, MYID;
	extern int WRITE_MODELFILES;
	extern char  MFILE[STRING_SIZE];	

	/* local variables */
	float rhov, muv, piv, vp, vs;
	float *pts, ts, tp, sumu, sumpi, ws;
	int i, j, l, ii, jj;
	char modfile[STRING_SIZE];	

	/*-----------------material property definition -------------------------*/	

	/* parameters for layer 1 */
	const float vp1=5100.0, vs1=2800.0, rho1=2100.0;


	/*-----------------------------------------------------------------------*/


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

			vp=vp1; vs=vs1; rhov=rho1; 


			muv=vs*vs*rhov/(1.0+sumu);
			piv=vp*vp*rhov/(1.0+sumpi);

			/* only the PE which belongs to the current global gridpoint
				  is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) &&
					(POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				taus[jj][ii]=ts;
				taup[jj][ii]=tp;
				u[jj][ii]=muv;
				rho[jj][ii]=rhov;
				pi[jj][ii]=piv;
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
}
