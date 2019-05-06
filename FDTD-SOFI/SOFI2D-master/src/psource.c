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
 *   generate P-wave source at source nodes
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float ** sxx, float ** syy,
		float **  srcpos_loc, float ** signals, int nsrc){

	extern float DH;
//	extern int RSG;
	extern int NT;
	int i, j, l;
	float amp=0;



	/* adding source wavelet to stress components 
	   (explosive source) at source points */

/*	if (RSG){
		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];

			//amp=signals[l][nt]/4.0; //unscaled explosive source
			//amp=(signals[l][nt])/(4.0*DH*DH); //scaled explosive source, seismic Moment = 1 Nm
			if(nt==1){amp=signals[l][nt+1]/(2.0*DH*DH);}
	                if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DH*DH);}
        	        if(nt==NT){amp=-signals[l][nt-1]/(2.0*DH*DH);}


			sxx[j][i]+=amp;
			sxx[j][i+1]+=amp;
			sxx[j+1][i]+=amp;
			sxx[j+1][i+1]+=amp;

			syy[j][i]+=amp;
			syy[j][i+1]+=amp;
			syy[j+1][i]+=amp;
			syy[j+1][i+1]+=amp;
		}

	}else{*/

		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];

			//amp=signals[l][nt]; //unscaled explosive source
			amp=(signals[l][nt])/(DH*DH); //scaled explosive source, seismic Moment = 1 Nm
			
			if(nt==1){amp=signals[l][nt+1]/(2.0*DH*DH);}
                        if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DH*DH);}
                        if(nt==NT){amp=-signals[l][nt-1]/(2.0*DH*DH);}

			
			sxx[j][i]+=amp;
			syy[j][i]+=amp;
		
	}
}
