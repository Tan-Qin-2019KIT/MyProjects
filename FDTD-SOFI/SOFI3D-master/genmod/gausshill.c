/*
 *   homogeneous space with a tunnel
 *   last update 05.02.99, T. Bohlen
 */

#include "fd.h"

void model(float  ***  rho, float ***  pi, float ***  u, 
float ***  taus, float ***  taup, float *  eta){

	/*--------------------------------------------------------------------------*/
	/* extern variables */
 
	extern float DT, DH, *FL, TS, TAU, REFREC[4];
	extern int NX, NY, NZ, NXG, NYG, NZG, POS[4], L, MYID;
	extern char  MFILE[STRING_SIZE];
	
	/* local variables */
	float muv, piv, ws;
	float Vp, Vs, Rho;
	float x, y, z, r, yhill;
	float *pts, ts, tp, sumu, sumpi;
	int i, j, k, l, ii, jj, kk;
	

	/* air */
	const float vp1=0.0, vs1=0.0, rho1=1.25, Qp1=10000.0, Qs1=10000.0;
/*	const float vp1=3000.0, vs1=1500.0, rho1=1200.0, Qp1=10000.0, Qs1=10000.0;*/

	/* hill */
	const float h=1000.0, a=1000.0, x0=REFREC[1], z0=REFREC[3], y0=REFREC[2], 
	vp2=3000.0, vs2=1500.0, rho2=1200.0;
	

	
	/*-----------------------------------------------------------------------*/


	/* vector for maxwellbodies */
	pts=vector(1,L);
	for (l=1;l<=L;l++) {
		pts[l]=1.0/(2.0*PI*FL[l]);
		eta[l]=DT/pts[l];
	}


	ws=2.0*PI*FL[1];



	/* loop over global grid */
	for (k=1;k<=NZG;k++){
		for (i=1;i<=NXG;i++){
		
 			for (j=1;j<=NYG;j++){
				
				x=(float)i*DH-x0;
				y=(float)j*DH;
				z=(float)k*DH-z0;


				r=sqrt(x*x+z*z);
				yhill=h*(1-exp(-(r*r)/(a*a)))+y0;
				if (y<yhill){
					Vp=vp1; Vs=vs1; Rho=rho1; tp=2.0/Qp1; ts=2.0/Qs1;}
				else {
					Vp=vp2; Vs=vs2; Rho=rho2; tp=TAU; ts=TAU;}
				
				


				sumu=0.0; 
				sumpi=0.0;
				for (l=1;l<=L;l++){
					sumu=sumu+((ws*ws*pts[l]*pts[l]*ts)/(1.0+ws*ws*pts[l]*pts[l]));
					sumpi=sumpi+((ws*ws*pts[l]*pts[l]*tp)/(1.0+ws*ws*pts[l]*pts[l]));
				}

				muv=Vs*Vs*Rho/(1.0+sumu);
				piv=Vp*Vp*Rho/(1.0+sumpi);

				/* only the PE which belongs to the current global gridpoint 
				  is saving model parameters in his local arrays */
				if ((POS[1]==((i-1)/NX)) && 
				    (POS[2]==((j-1)/NY)) && 
				    (POS[3]==((k-1)/NZ))){
					ii=i-POS[1]*NX;
					jj=j-POS[2]*NY;
					kk=k-POS[3]*NZ;

					taus[jj][ii][kk]=ts;
					taup[jj][ii][kk]=tp;
					u[jj][ii][kk]=muv;
					rho[jj][ii][kk]=Rho;
					pi[jj][ii][kk]=piv;
				}
			}
		}
	}	

 

/*	writemod(MFILE,rho,3);

	MPI_Barrier(MPI_COMM_WORLD);

	if (MYID==0) mergemod(MFILE,3);
*/
	free_vector(pts,1,L);
	
}



