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
 *   updating particle velocities at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of FDORDER accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"



void update_v(int nx1, int nx2, int ny1, int ny2, int nt,
		float **  vx, float ** vy, float ** sxx, float ** syy,
		float ** sxy, float **rho, float  **rip, float **rjp,
		float **  srcpos_loc, float ** signals, int nsrc, float ** absorb_coeff,
		float *hc){

	int i, j,l;
	//int m, fdoh;
	float amp, dtdh, dtdh2;
	//float vxtmp, vytmp;
	float ripjm;
	float sxx_xs, sxx_ys, sxy_xs, sxy_ys, syy_xs, syy_ys, sxx_x, sxy_x, sxy_y, syy_y;
	float azi_rad;

	extern float DT, DH, DH;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;
	extern int MYID, SOURCE_TYPE, CHECKPTREAD, FDORDER;
	//extern int NPROCX, NPROCY, POS[3];
	extern FILE *FP;
	extern int RSG;

	//fdoh = FDORDER/2;
	dtdh = DT/DH;
	dtdh2=dtdh/2.0;


	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_v (printed by PE %d):\n",MYID);
		fprintf(FP," Updating particle velocities ...");
	}


	/* ------------------------------------------------------------
	 * Important!
	 * rip and rjp are reciprocal values of averaged densities
	 * ------------------------------------------------------------ */



	if (!CHECKPTREAD)
		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
			azi_rad=srcpos_loc[7][l]*PI/180;

			//amp=signals[l][nt]; // unscaled force amplitude
			amp=(DT*signals[l][nt])/(DH*DH);// scaled force amplitude with F= 1N

			//fprintf(FP," amp at timestep nt %i = %5.5e with DH=%5.2f  DT=%5.8f\n",nt,amp,DH,DT);

			SOURCE_TYPE=(int)srcpos_loc[8][l];

			switch (SOURCE_TYPE){
			case 2 : /* single force in x */


				vx[j][i]  +=  rip[j][i]*amp;

				/* previous implementation of body forces as seismic sources.
				 * Implementation according to Coutant et al., BSSA, Vol. 85, No 5, 1507-1512.
				 * The stress tensor components sxx and syy are incremented prior to
				 * particle velocity update. Thereby the body force (both directions)
				 * are located at full grid point (i,j) (same position as pressure source).
				 * as a consequence, source signals are added [and weighted] at multiple grid points.
				 * This implementation works but not quite physical when considering e.g. a force of 1 N
				 * and aiming to gain the particle velocity strictly according to that force.
				 * Therefore it has been commented */

				/*for (m=1; m<=fdoh; m++) {

					vx[j][i+m-1]  +=  hc[m]*rip[j][i]*amp;
					vx[j][i-m]    +=  hc[m]*rip[j][i-1]*amp;

				}*/


				break;
			case 3 : /* single force in y */

				vy[j][i]  +=  rjp[j][i]*amp;

				/*for (m=1; m<=fdoh; m++) {

					vy[j+m-1][i]  +=  hc[m]*rjp[j][i]*amp;
					vy[j-m][i]    +=  hc[m]*rjp[j][i-1]*amp;

				}*/

				break;
			case 4 : /* custom force */

				vx[j][i]  +=  sin(azi_rad)*(rip[j][i]*amp);
				vy[j][i]  +=  cos(azi_rad)*(rjp[j][i]*amp);

				/*for (m=1; m<=fdoh; m++) {

					vx[j][i+m-1]  +=  sin(azi_rad)*(hc[m]*rip[j][i]*amp);
					vx[j][i-m]    +=  sin(azi_rad)*(hc[m]*rip[j][i-1]*amp);

					vy[j+m-1][i]  +=  cos(azi_rad)*(hc[m]*rjp[j][i]*amp);
					vy[j-m][i]    +=  cos(azi_rad)*(hc[m]*rjp[j][i-1]*amp);

				}*/

				break;

			}
		}



	if (RSG) {/* rotated staggered grid of second order accuracy (time and space) */
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

				/* interpolate density between gridpoints*/
				ripjm=0.25*(rho[j][i]+rho[j][i-1]+rho[j-1][i]+rho[j-1][i-1]); 	

				/* 2th order, in the rotated coordinate system: */
				sxx_xs=sxx[j-1][i]-sxx[j][i-1];
				sxx_ys=sxx[j][i]-sxx[j-1][i-1];
				sxy_xs=sxy[j-1][i]-sxy[j][i-1];
				sxy_ys=sxy[j][i]-sxy[j-1][i-1];
				syy_xs=syy[j-1][i]-syy[j][i-1];
				syy_ys=syy[j][i]-syy[j-1][i-1];

				/* in the cartesian coordinate system: */
				sxx_x=sxx_ys+sxx_xs;
				sxy_x=sxy_ys+sxy_xs;
				sxy_y=sxy_ys-sxy_xs;
				syy_y=syy_ys-syy_xs;			

				vx[j][i]+= dtdh2*(sxx_x+sxy_y)/ripjm;
				vy[j][i]+= dtdh2*(sxy_x+syy_y)/ripjm;
			}
		}
	}else{


	switch (FDORDER){
	case 2:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

				/* updating the x-component of the velocity (vx) */
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i])	+ hc[1]*(sxy[j][i]-sxy[j-1][i])	)*dtdh*rip[j][i];

				/* updating the y-component of the velocity (vy) */
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i])	+ hc[1]*(sxy[j][i]-sxy[j][i-1])	)*dtdh*rjp[j][i];

			}
		}
		break;

	case 4:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
						+ hc[2]*(sxx[j][i+2]-sxx[j][i-1])
						+ hc[1]*(sxy[j][i]  -sxy[j-1][i])
						+ hc[2]*(sxy[j+1][i]-sxy[j-2][i])
				)*dtdh*rip[j][i];

				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
						+ hc[2]*(syy[j+2][i]-syy[j-1][i])
						+ hc[1]*(sxy[j][i]  -sxy[j][i-1])
						+ hc[2]*(sxy[j][i+1]-sxy[j][i-2])
				)*dtdh*rjp[j][i];
			}
		}
		break;

	case 6:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
						+ hc[2]*(sxx[j][i+2]-sxx[j][i-1])
						+ hc[3]*(sxx[j][i+3]-sxx[j][i-2])
						+ hc[1]*(sxy[j][i]  -sxy[j-1][i])
						+ hc[2]*(sxy[j+1][i]-sxy[j-2][i])
						+ hc[3]*(sxy[j+2][i]-sxy[j-3][i])
				)*dtdh*rip[j][i];

				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
						+ hc[2]*(syy[j+2][i]-syy[j-1][i])
						+ hc[3]*(syy[j+3][i]-syy[j-2][i])
						+ hc[1]*(sxy[j][i]  -sxy[j][i-1])
						+ hc[2]*(sxy[j][i+1]-sxy[j][i-2])
						+ hc[3]*(sxy[j][i+2]-sxy[j][i-3])
				)*dtdh*rjp[j][i];
			}
		}
		break;

	case 8:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			  
			  
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
						+ hc[2]*(sxx[j][i+2]-sxx[j][i-1])
						+ hc[3]*(sxx[j][i+3]-sxx[j][i-2])
						+ hc[4]*(sxx[j][i+4]-sxx[j][i-3])
						+ hc[1]*(sxy[j][i]  -sxy[j-1][i])
						+ hc[2]*(sxy[j+1][i]-sxy[j-2][i])
						+ hc[3]*(sxy[j+2][i]-sxy[j-3][i])
						+ hc[4]*(sxy[j+3][i]-sxy[j-4][i])
				)*dtdh*rip[j][i];

				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
						+ hc[2]*(syy[j+2][i]-syy[j-1][i])
						+ hc[3]*(syy[j+3][i]-syy[j-2][i])
						+ hc[4]*(syy[j+4][i]-syy[j-3][i])
						+ hc[1]*(sxy[j][i]  -sxy[j][i-1])
						+ hc[2]*(sxy[j][i+1]-sxy[j][i-2])
						+ hc[3]*(sxy[j][i+2]-sxy[j][i-3])
						+ hc[4]*(sxy[j][i+3]-sxy[j][i-4])
				)*dtdh*rjp[j][i];
			}
		}
		break;

	case 10:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
						+ hc[2]*(sxx[j][i+2]-sxx[j][i-1])
						+ hc[3]*(sxx[j][i+3]-sxx[j][i-2])
						+ hc[4]*(sxx[j][i+4]-sxx[j][i-3])
						+ hc[5]*(sxx[j][i+5]-sxx[j][i-4])
						+ hc[1]*(sxy[j][i]  -sxy[j-1][i])
						+ hc[2]*(sxy[j+1][i]-sxy[j-2][i])
						+ hc[3]*(sxy[j+2][i]-sxy[j-3][i])
						+ hc[4]*(sxy[j+3][i]-sxy[j-4][i])
						+ hc[5]*(sxy[j+4][i]-sxy[j-5][i])
				)*dtdh*rip[j][i];

				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
						+ hc[2]*(syy[j+2][i]-syy[j-1][i])
						+ hc[3]*(syy[j+3][i]-syy[j-2][i])
						+ hc[4]*(syy[j+4][i]-syy[j-3][i])
						+ hc[5]*(syy[j+5][i]-syy[j-4][i])
						+ hc[1]*(sxy[j][i]  -sxy[j][i-1])
						+ hc[2]*(sxy[j][i+1]-sxy[j][i-2])
						+ hc[3]*(sxy[j][i+2]-sxy[j][i-3])
						+ hc[4]*(sxy[j][i+3]-sxy[j][i-4])
						+ hc[5]*(sxy[j][i+4]-sxy[j][i-5])
				)*dtdh*rjp[j][i];
			}
		}
		break;

	case 12:
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i]  )
						+ hc[2]*(sxx[j][i+2]-sxx[j][i-1])
						+ hc[3]*(sxx[j][i+3]-sxx[j][i-2])
						+ hc[4]*(sxx[j][i+4]-sxx[j][i-3])
						+ hc[5]*(sxx[j][i+5]-sxx[j][i-4])
						+ hc[6]*(sxx[j][i+6]-sxx[j][i-5])
						+ hc[1]*(sxy[j][i]  -sxy[j-1][i])
						+ hc[2]*(sxy[j+1][i]-sxy[j-2][i])
						+ hc[3]*(sxy[j+2][i]-sxy[j-3][i])
						+ hc[4]*(sxy[j+3][i]-sxy[j-4][i])
						+ hc[5]*(sxy[j+4][i]-sxy[j-5][i])
						+ hc[6]*(sxy[j+5][i]-sxy[j-6][i])
				)*dtdh*rip[j][i];

				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i]  )
						+ hc[2]*(syy[j+2][i]-syy[j-1][i])
						+ hc[3]*(syy[j+3][i]-syy[j-2][i])
						+ hc[4]*(syy[j+4][i]-syy[j-3][i])
						+ hc[5]*(syy[j+5][i]-syy[j-4][i])
						+ hc[6]*(syy[j+6][i]-syy[j-5][i])
						+ hc[1]*(sxy[j][i]  -sxy[j][i-1])
						+ hc[2]*(sxy[j][i+1]-sxy[j][i-2])
						+ hc[3]*(sxy[j][i+2]-sxy[j][i-3])
						+ hc[4]*(sxy[j][i+3]-sxy[j][i-4])
						+ hc[5]*(sxy[j][i+4]-sxy[j][i-5])
						+ hc[6]*(sxy[j][i+5]-sxy[j][i-6])
				)*dtdh*rjp[j][i];
			}
		}
		break;

	default: // 2nd order
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

				/* vxtmp = 0;
				vytmp = 0;
				for (m=1; m<=fdoh; m++) {
					vxtmp +=   hc[m]*( sxx[j][i+m]   - sxx[j][i-m+1] )
														 + hc[m]*( sxy[j+m-1][i] - sxy[j-m][i]   );

					vytmp +=   hc[m]*( syy[j+m][i]   - syy[j-m+1][i] )
														 + hc[m]*( sxy[j][i+m-1] - sxy[j][i-m]   );
					vx[j][i] += vxtmp*dtdh*rip[j][i];
					vy[j][i] += vytmp*dtdh*rjp[j][i];

				}*/

				/* updating the x-component of the velocity (vx) */
				vx[j][i]+= (  hc[1]*(sxx[j][i+1]-sxx[j][i])	+ hc[1]*(sxy[j][i]-sxy[j-1][i]))*dtdh*rip[j][i];

				/* updating the y-component of the velocity (vy) */
				vy[j][i]+= (  hc[1]*(syy[j+1][i]-syy[j][i])	+ hc[1]*(sxy[j][i]-sxy[j][i-1]))*dtdh*rjp[j][i];

			}
		}
		break;

	} /* end of switch(FDORDER) */
	}



	/*if ((FW) && (ABS_TYPE == 2))
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){
			           vx[j][i]*=absorb_coeff[j][i];
				vy[j][i]*=absorb_coeff[j][i];
			       
				sxy[j][i]*=absorb_coeff[j][i];
				sxx[j][i]*=absorb_coeff[j][i];
				syy[j][i]*=absorb_coeff[j][i];
			}
		} */

	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.3f s).\n",time2-time1);
	}
}
