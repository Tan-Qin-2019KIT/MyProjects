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
 *   store amplitudes (particle velocities or pressure or curl and div) 
 *   at receiver positions in arrays
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo_ssg(int lsamp, int ntr, int **recpos, float **sectionvx, 
float **sectionvy, float **sectionp, float **sectioncurl, float **sectiondiv,
float **vx, float **vy, float **sxx, float **syy, float **pi, float **u, float *hc){ 
		
	extern int NDT, SEISMO, FDORDER;	
	extern float DH;
	int i,j, itr, ins, nxrec, nyrec, m, fdoh;
	float dhi, vxx, vyy, vxy, vyx;
	//float dh24;


	dhi = 1.0/DH;
	//dh24=1.0/(DH*24.0);
	fdoh = FDORDER/2;


	ins=lsamp/NDT;
	for (itr=1;itr<=ntr;itr++){
		nxrec=recpos[1][itr];
		nyrec=recpos[2][itr];
		switch (SEISMO){
		case 1 : /* particle velocities */
			sectionvx[itr][ins]=vx[nyrec][nxrec];
			sectionvy[itr][ins]=vy[nyrec][nxrec];
			break;
		
		case 2 : /* pressure */
			i=nxrec; j=nyrec;
			//sectionp[itr][ins]=-sxx[nyrec][nxrec]-syy[nyrec][nxrec]; // unscaled amplitude
			sectionp[itr][ins]=((3.0*pi[j][i]-4.0*u[j][i])/(2.0*pi[j][i]-2.0*u[j][i]))*(-sxx[nyrec][nxrec]-syy[nyrec][nxrec])/3; // true amplitude
			break;
		
		case 3 : /* curl +div */
			i=nxrec; j=nyrec;
			
			vxx = 0;
			vyy = 0;
			vyx = 0;
			vxy = 0;
			for (m=1; m<=fdoh; m++) {
				vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
				vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
			}
			vxx *= dhi;
			vyy *= dhi;
			vyx *= dhi;
			vxy *= dhi;
			
			sectiondiv[itr][ins]=(vxx+vyy)*sqrt(pi[j][i]);
			sectioncurl[itr][ins]=(vxy-vyx)*sqrt(u[j][i]);
			break;
		
		case 4 : /* all */
			i=nxrec; j=nyrec;

			vxx = 0;
			vyy = 0;
			vyx = 0;
			vxy = 0;
			for (m=1; m<=fdoh; m++) {
				vxx += hc[m]*(vx[j][i+m-1] -vx[j][i-m]  );
				vyy += hc[m]*(vy[j+m-1][i] -vy[j-m][i]  );
				vyx += hc[m]*(vy[j][i+m]   -vy[j][i-m+1]);
				vxy += hc[m]*(vx[j+m][i]   -vx[j-m+1][i]);
			}
			vxx *= dhi;
			vyy *= dhi;
			vyx *= dhi;
			vxy *= dhi;

			sectiondiv[itr][ins]=(vxx+vyy)*sqrt(pi[j][i]);
			sectioncurl[itr][ins]=(vxy-vyx)*sqrt(u[j][i]);
			sectionvx[itr][ins]=vx[nyrec][nxrec];
			sectionvy[itr][ins]=vy[nyrec][nxrec];			
			//sectionp[itr][ins]=-sxx[nyrec][nxrec]-syy[nyrec][nxrec]; //unscaled amplitude
			sectionp[itr][ins]=((3.0*pi[j][i]-4.0*u[j][i])/(2.0*pi[j][i]-2.0*u[j][i]))*(-sxx[nyrec][nxrec]-syy[nyrec][nxrec])/3;
			break;

		}

	}
}
