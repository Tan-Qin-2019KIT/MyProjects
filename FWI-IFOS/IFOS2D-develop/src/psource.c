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
 *   generate P-wave source at source nodes
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float ** sxx, float ** syy, float ** sp,
float **  srcpos_loc, float ** signals, int nsrc, int sw){


        extern int MYID, NT, ACOUSTIC;
	extern float DH;
	int i, j, l;
	float amp;



	/* adding source wavelet to stress components 
	   (explosive source) at source points */


	/* Forward Modelling (sw==0) */
	if(sw==0){
	for (l=1;l<=nsrc;l++) {
		i=(int)srcpos_loc[1][l];
		j=(int)srcpos_loc[2][l];
		
		if(nt==1){amp=signals[l][nt+1]/(2.0*DH*DH);}
		if((nt>1)&&(nt<NT)){amp=(signals[l][nt+1]-signals[l][nt-1])/(2.0*DH*DH);}
		if(nt==NT){amp=-signals[l][nt-1]/(2.0*DH*DH);}
	
// 		amp = signals[l][nt]; //unscaled explosive source
		
		if(!ACOUSTIC){
			sxx[j][i]+=0.5*amp;
			syy[j][i]+=0.5*amp;
		}else
			sp[j][i]+=amp;
	}}

	/* Backpropagation (sw==1) */
	if(sw==1){
		for (l=1;l<=nsrc;l++) {
			i=(int)srcpos_loc[1][l];
			j=(int)srcpos_loc[2][l];
					
			amp = signals[l][nt];
			
			if(!ACOUSTIC){
				sxx[j][i]+=0.5*amp;
				syy[j][i]+=0.5*amp;
			}else
				sp[j][i]+=amp;
		}
	}
}
