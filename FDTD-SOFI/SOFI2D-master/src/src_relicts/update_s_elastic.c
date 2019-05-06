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
 *   updating stress components at gridpoints [nx1...nx2][ny1...ny2]
 *   by a staggered grid finite difference scheme of arbitrary (FDORDER) order accuracy in space
 *   and second order accuracy in time
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void update_s_elastic(int nx1, int nx2, int ny1, int ny2, int nt,
		float **  vx, float **   vy, float **   sxx, float **   syy,
		float **   sxy, float ** pi, float ** u, float ** uipjp, float ** absorb_coeff,
		float *hc){


	int i,j, m, fdoh;
	float fipjp, f, g;
	float vxxs, vxys, vyxs, vyys;
	float  vxx, vyy, vxy, vyx;
	float  dhi, dhi2;
	extern float DT, DH;
	extern int MYID, FDORDER, RSG;
	extern FILE *FP;
	extern int OUTNTIMESTEPINFO;
	double time1=0.0, time2=0.0;



	dhi = 1.0/DH;
	dhi2=dhi/2.0;
	fdoh = FDORDER/2;

	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
		time1=MPI_Wtime();
		fprintf(FP,"\n **Message from update_s (printed by PE %d):\n",MYID);
		fprintf(FP," Updating stress components ...");
	}	



	if (RSG)
		for (j=ny1;j<=ny2;j++){
			for (i=nx1;i<=nx2;i++){

				/* spatial derivatives of the components of particle velocities */

				/* in the rotated coordinate system: */
				/* second order */
				vxxs=vx[j][i+1]-vx[j+1][i];
				vxys=vx[j+1][i+1]-vx[j][i];
				vyxs=vy[j][i+1]-vy[j+1][i];
				vyys=vy[j+1][i+1]-vy[j][i];


				/* in the cartesian coordinate system: */
				vxx=(vxys+vxxs)*dhi2;
				vyx=(vyys+vyxs)*dhi2;
				vxy=(vxys-vxxs)*dhi2;
				vyy=(vyys-vyxs)*dhi2;




				/* updating components of the stress tensor, partially */
				f=u[j][i]*DT;
				g=pi[j][i]*DT;

				sxy[j][i]+=(f*(vxy+vyx));
				sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
				syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);

			}
		}/* end of if RSG=1 */
	else {

		switch (FDORDER){

		case 2:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
				
				                
				  /* Compute values for shearmodulus u[j][i],
				    P-wave modulus pi[j][i],
				    tau for S-waves and P-waves taus[j][i], 
				    taup[j][i] at staggered grid points: */

					/* spatial derivatives of the components of the velocities */
					/* using Holberg coefficients */
					vxx = hc[1]*(vx[j][i]  -vx[j][i-1])*dhi;
					vyy = hc[1]*(vy[j][i]  -vy[j-1][i])*dhi;
					vyx = hc[1]*(vy[j][i+1]-vy[j][i])*dhi;
					vxy = hc[1]*(vx[j+1][i]-vx[j][i])*dhi;

					/* updating components of the stress tensor, partially */
					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				}
			}
			break;

		case 4:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
							+ hc[2]*(vx[j][i+1]-vx[j][i-2])
					)*dhi;
					vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
							+ hc[2]*(vy[j+1][i]-vy[j-2][i])
					)*dhi;
					vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
							+ hc[2]*(vy[j][i+2]-vy[j][i-1])
					)*dhi;
					vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
							+ hc[2]*(vx[j+2][i]-vx[j-1][i])
					)*dhi;

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				}
			}
			break;

		case 6:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
							+ hc[2]*(vx[j][i+1]-vx[j][i-2])
							+ hc[3]*(vx[j][i+2]-vx[j][i-3])
					)*dhi;
					vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
							+ hc[2]*(vy[j+1][i]-vy[j-2][i])
							+ hc[3]*(vy[j+2][i]-vy[j-3][i])
					)*dhi;
					vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
							+ hc[2]*(vy[j][i+2]-vy[j][i-1])
							+ hc[3]*(vy[j][i+3]-vy[j][i-2])
					)*dhi;
					vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
							+ hc[2]*(vx[j+2][i]-vx[j-1][i])
							+ hc[3]*(vx[j+3][i]-vx[j-2][i])
					)*dhi;

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				}
			}
			break;

		case 8:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
							+ hc[2]*(vx[j][i+1]-vx[j][i-2])
							+ hc[3]*(vx[j][i+2]-vx[j][i-3])
							+ hc[4]*(vx[j][i+3]-vx[j][i-4])
					)*dhi;
					vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
							+ hc[2]*(vy[j+1][i]-vy[j-2][i])
							+ hc[3]*(vy[j+2][i]-vy[j-3][i])
							+ hc[4]*(vy[j+3][i]-vy[j-4][i])
					)*dhi;
					vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
							+ hc[2]*(vy[j][i+2]-vy[j][i-1])
							+ hc[3]*(vy[j][i+3]-vy[j][i-2])
							+ hc[4]*(vy[j][i+4]-vy[j][i-3])
					)*dhi;
					vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
							+ hc[2]*(vx[j+2][i]-vx[j-1][i])
							+ hc[3]*(vx[j+3][i]-vx[j-2][i])
							+ hc[4]*(vx[j+4][i]-vx[j-3][i])
					)*dhi;

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				}
			}
			break;

		case 10:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
							+ hc[2]*(vx[j][i+1]-vx[j][i-2])
							+ hc[3]*(vx[j][i+2]-vx[j][i-3])
							+ hc[4]*(vx[j][i+3]-vx[j][i-4])
							+ hc[5]*(vx[j][i+4]-vx[j][i-5])
					)*dhi;
					vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
							+ hc[2]*(vy[j+1][i]-vy[j-2][i])
							+ hc[3]*(vy[j+2][i]-vy[j-3][i])
							+ hc[4]*(vy[j+3][i]-vy[j-4][i])
							+ hc[5]*(vy[j+4][i]-vy[j-5][i])
					)*dhi;
					vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
							+ hc[2]*(vy[j][i+2]-vy[j][i-1])
							+ hc[3]*(vy[j][i+3]-vy[j][i-2])
							+ hc[4]*(vy[j][i+4]-vy[j][i-3])
							+ hc[5]*(vy[j][i+5]-vy[j][i-4])
					)*dhi;
					vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
							+ hc[2]*(vx[j+2][i]-vx[j-1][i])
							+ hc[3]*(vx[j+3][i]-vx[j-2][i])
							+ hc[4]*(vx[j+4][i]-vx[j-3][i])
							+ hc[5]*(vx[j+5][i]-vx[j-4][i])
					)*dhi;

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);

					/*
				piv = pi[j][i]*(vxx+vyy);
				sxy[j][i]+=(uipjp[j][i]*(vxy+vyx));
				sxx[j][i]+=piv-(2.0*u[j][i]*vyy);
				syy[j][i]+=piv-(2.0*u[j][i]*vxx);
					 */
				}
			}
			break;

		case 12:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = (  hc[1]*(vx[j][i]  -vx[j][i-1])
							+ hc[2]*(vx[j][i+1]-vx[j][i-2])
							+ hc[3]*(vx[j][i+2]-vx[j][i-3])
							+ hc[4]*(vx[j][i+3]-vx[j][i-4])
							+ hc[5]*(vx[j][i+4]-vx[j][i-5])
							+ hc[6]*(vx[j][i+5]-vx[j][i-6])
					)*dhi;
					vyy = (  hc[1]*(vy[j][i]  -vy[j-1][i])
							+ hc[2]*(vy[j+1][i]-vy[j-2][i])
							+ hc[3]*(vy[j+2][i]-vy[j-3][i])
							+ hc[4]*(vy[j+3][i]-vy[j-4][i])
							+ hc[5]*(vy[j+4][i]-vy[j-5][i])
							+ hc[6]*(vy[j+5][i]-vy[j-6][i])
					)*dhi;
					vyx = (  hc[1]*(vy[j][i+1]-vy[j][i])
							+ hc[2]*(vy[j][i+2]-vy[j][i-1])
							+ hc[3]*(vy[j][i+3]-vy[j][i-2])
							+ hc[4]*(vy[j][i+4]-vy[j][i-3])
							+ hc[5]*(vy[j][i+5]-vy[j][i-4])
							+ hc[6]*(vy[j][i+6]-vy[j][i-5])
					)*dhi;
					vxy = (  hc[1]*(vx[j+1][i]-vx[j][i])
							+ hc[2]*(vx[j+2][i]-vx[j-1][i])
							+ hc[3]*(vx[j+3][i]-vx[j-2][i])
							+ hc[4]*(vx[j+4][i]-vx[j-3][i])
							+ hc[5]*(vx[j+5][i]-vx[j-4][i])
							+ hc[6]*(vx[j+6][i]-vx[j-5][i])
					)*dhi;

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx));
					sxx[j][i]+=(g*(vxx+vyy))-(2.0*f*vyy);
					syy[j][i]+=(g*(vxx+vyy))-(2.0*f*vxx);
				}
			}
			break;

		default:
			for (j=ny1;j<=ny2;j++){
				for (i=nx1;i<=nx2;i++){
					vxx = 0.0;
					vyy = 0.0;
					vyx = 0.0;
					vxy = 0.0;
					for (m=1; m<=fdoh; m++) {
						vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
						vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
						vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
						vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
					}

					fipjp=uipjp[j][i]*DT;
					f=u[j][i]*DT;
					g=pi[j][i]*DT;

					sxy[j][i]+=(fipjp*(vxy+vyx))*dhi;
					sxx[j][i]+=((g*(vxx+vyy))-(2.0*f*vyy))*dhi;
					syy[j][i]+=((g*(vxx+vyy))-(2.0*f*vxx))*dhi;
				}
			}
			break;

		} /* end of switch(FDORDER) */
	}/* end of if RSG!=1 */


	if ((MYID==0) && ((nt+(OUTNTIMESTEPINFO-1))%OUTNTIMESTEPINFO)==0) {
		time2=MPI_Wtime();
		fprintf(FP," finished (real time: %4.3f s).\n",time2-time1);
	}
}
