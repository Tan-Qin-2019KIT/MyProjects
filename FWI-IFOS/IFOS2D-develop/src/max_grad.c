/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   calculate test step length for material parameter update
 *   
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void max_grad(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float **  pi, float **  u){


	/*--------------------------------------------------------------------------*/
	/* extern variables */
	extern float DH, DT;
	extern float EPSILON, EPSILON_u, EPSILON_rho;
	extern int NX, NY, NXG, NYG,  POS[3], MYID;


	/* local variables */

	float  Zp, Zs;
	float  pimax, rhomax, umax, gradmax, gradmax_rho, gradmax_u, epsilon1, pimaxr, gradmaxr, gradmaxr_u, umaxr, gradmaxr_rho, rhomaxr;
	int i, j;


        /* invert for Zp and Zs */
	/* ------------------------------------------------------------------------------------ */

	
	/* find maximum of Zp and gradient waveconv */
	pimax = 0.0;
	gradmax = 0.0;

    for (j=1;j<=NY;j++){
	    for (i=1;i<=NX;i++){
		
		Zp = sqrt((pi[j][i] + 2.0 * u[j][i])*rho[j][i]);
		
		
		if(Zp>pimax){pimax=Zp;}
		
		if((i*j == 1) || (gradmax == 0.0)) {
				gradmax = fabs(waveconv[j][i]);		
		} else {		
		   if(fabs(waveconv[j][i]) > gradmax){
		      gradmax = fabs(waveconv[j][i]);
		   }
		                               
		}
	}}
	
	/* find maximum of Zs and gradient waveconv_u */
	
	umax = 0.0;
	gradmax_u = 0.0;
	
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
		
		Zs = sqrt(u[j][i]*rho[j][i]);
		
		if(Zs>umax){umax=Zs;}
		
		if((i*j == 1) || (gradmax_u == 0.0)) {
		        gradmax_u = fabs(waveconv_u[j][i]);
		} else {
		if(fabs(waveconv_u[j][i]) > gradmax_u)
		{
		gradmax_u = fabs(waveconv_u[j][i]);
					}
				}
		                               
		}
	}
	
	
	
	
	/* find maximum of rho and gradient waveconv_rho */
	rhomax = 0.0;
	gradmax_rho = 0.0;
	
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
		
		if(rho[j][i]>rhomax){rhomax=rho[j][i];}
		
		if((i*j == 1) || (gradmax_rho == 0.0)) {
		        gradmax_rho = fabs(waveconv_rho[j][i]);
		} else {
		if(fabs(waveconv_rho[j][i]) > gradmax_rho)
		{
		gradmax_rho = fabs(waveconv_rho[j][i]);
					}
				}
		                               
		}
	}
		
	/* calculate scaling factor for the gradients */
        /* --------------------------------------------- */
	
	/* parameter 1 */
	   MPI_Allreduce(&pimax,  &pimaxr,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax,&gradmaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	
	   EPSILON = (pimaxr/gradmaxr);
	   epsilon1=EPSILON;
	   MPI_Allreduce(&EPSILON,&epsilon1,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	   if (MYID==0)  EPSILON=epsilon1;

	   exchange_par();
	

	
	/* parameter 2 */
	   MPI_Allreduce(&umax,&umaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax_u,&gradmaxr_u,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
	   EPSILON_u = (umaxr/gradmaxr_u);

           epsilon1=EPSILON_u;
	   MPI_Allreduce(&EPSILON_u,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	   if (MYID==0)  EPSILON_u=epsilon1;

	   exchange_par();		
	
	/* parameter 3 */

	   MPI_Allreduce(&rhomax,&rhomaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(&gradmax_rho,&gradmaxr_rho,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	
    	   EPSILON_rho = (rhomaxr/gradmaxr_rho);
           epsilon1=EPSILON_rho;
	   MPI_Allreduce(&EPSILON_rho,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
	
	   if (MYID==0)  EPSILON_rho=epsilon1;

	   exchange_par();		

        /* apply scaling factors */
	/* --------------------- */
	
    for (j=1;j<=NY;j++){
        for (i=1;i<=NX;i++){
		
		waveconv[j][i] = EPSILON * waveconv[j][i];
		waveconv_u[j][i] = EPSILON_u * waveconv_u[j][i];
		waveconv_rho[j][i] = EPSILON_rho * waveconv_rho[j][i];
		                              
		}
	}	   
	   
		
		
}

