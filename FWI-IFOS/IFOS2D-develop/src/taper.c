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
 *   Taper function now adapted for use in stf.c, only for single traces!
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  taper(float *section, int ns, float fc){

	/* declaration of extern variables */
	extern int MYID;
	extern float DT;
	extern FILE *FP;
	
	/* declaration of local variables */
	int i,j, h, taperlength, taperduration;
	int tracl1;
	float a;
	float damping, amp;
	float *window=NULL, *amp1=NULL;
	
	window = vector(1,ns);
        amp1 = vector(1,ns);
	
	taperlength=(int)(ceil(2.0/fc/DT));
	taperduration=2*taperlength;
	
	/* "Cerjan"-Window */
        damping=99.9;
        amp=1.0-damping/100.0;
	        a=sqrt(-log(amp)/((taperlength-1)*(taperlength-1)));
        
	for (i=1;i<=ns;i++){
		window[i]=1.0;
		amp1[i]=0.0;
	}
	
	if (MYID==0){
		fprintf(FP,"\n fc: %f\n",fc);
		fprintf(FP,"\n taperlength: %d\n",taperlength);
	}
	
	for (i=1;i<=taperlength;i++){
		amp1[i]=exp(-(a*a*(taperlength-i)*(taperlength-i)));
	}
	
// 	/* Taper at the beginning of the window*/
// 	for (i=1;i<=taperlength;i++){
// 		window[i]=amp1[i];
// 	}
	
	h=1;
	for (i=taperduration;i>=(taperduration-taperlength+3);i--){
		window[i]=amp1[h];
		h++;
	}
	
	for (i=taperduration;i<=ns;i++){
		window[i]=amp1[i];
	}
	
	for(j=1;j<=ns;j++){
		section[j]*=window[j];
	}
	
	free_vector(window,1,ns);
	free_vector(amp1,1,ns);
}
